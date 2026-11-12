"""Method agreement studies over the existing measurement ledger.

A study fixes:

* the two methods (instrument/channel pairs, reference and candidate);
* the specimen tag used to identify the same physical specimen;
* the replicate strategy (first, last or mean of repeated readings);
* the valid time window and the maximum within-pair time separation;
* the pre-determined acceptable difference and the claimed measurement range.

Pairing is key based, never positional: observations are matched through the
specimen key plus time proximity, every observation is consumed at most once,
and every exclusion is reported with an explicit reason.  Repeated readings of
one specimen are aggregated before statistics so a specimen is a single
independent unit.
"""

from dataclasses import dataclass, field
from decimal import localcontext

from ..errors import Conflict, InputError, Missing, StateError
from ..ledger.audit import append
from ..observations import Measurement
from ..values import (
    canonical,
    decimal_text,
    digest,
    identifier,
    instant,
    integer,
    number,
    require_object,
    strict_json,
    time_text,
)
from .statistics import (
    MAX_BOOTSTRAP,
    MAX_POINTS,
    Point,
    ba_bootstrap_statistic,
    bland_altman,
    cluster_bootstrap,
    deming,
    transform,
)

REPLICATE_STRATEGIES = ("first", "last", "mean")
SCALES = ("absolute", "percent", "log")
MIN_BOOTSTRAP = 10


@dataclass(frozen=True)
class StudySpec:
    id: str
    reference_instrument: str
    reference_channel: str
    candidate_instrument: str
    candidate_channel: str
    start: object
    end: object
    specimen_tag: str
    replicate_strategy: str
    pairing_seconds: int
    acceptable_difference: object
    range_low: object
    range_high: object
    scale: str
    percent_fraction: object
    bootstrap_replicates: int
    bootstrap_seed: int
    variance_ratio: object
    target_unit: str

    @classmethod
    def parse(cls, data):
        fields = {
            "id",
            "reference",
            "candidate",
            "start",
            "end",
            "specimen_tag",
            "replicate_strategy",
            "pairing_seconds",
            "acceptable_difference",
            "measurement_range",
            "scale",
            "percent_fraction",
            "bootstrap_replicates",
            "bootstrap_seed",
            "variance_ratio",
            "target_unit",
        }
        require_object(data, fields, {
            "id",
            "reference",
            "candidate",
            "start",
            "end",
            "specimen_tag",
            "replicate_strategy",
            "pairing_seconds",
            "acceptable_difference",
            "measurement_range",
        })
        reference = _method(data["reference"], "reference")
        candidate = _method(data["candidate"], "candidate")
        if reference == candidate:
            raise InputError("reference and candidate must be distinct methods")
        start, end = instant(data["start"]), instant(data["end"])
        if end <= start:
            raise InputError("study interval must increase")
        specimen_tag = identifier(data["specimen_tag"], "specimen_tag")
        strategy = data["replicate_strategy"]
        if strategy not in REPLICATE_STRATEGIES:
            raise InputError("unsupported replicate strategy", field="replicate_strategy")
        pairing_seconds = integer(
            data["pairing_seconds"], "pairing_seconds", 0, 86400
        )
        acceptable = number(
            data["acceptable_difference"], "acceptable_difference", positive=True
        )
        rng = data["measurement_range"]
        require_object(rng, {"low", "high"}, {"low", "high"})
        range_low, range_high = number(rng["low"], "range.low"), number(
            rng["high"], "range.high"
        )
        if range_high <= range_low:
            raise InputError("measurement range must increase")
        scale = data.get("scale", "absolute")
        if scale not in SCALES:
            raise InputError("unsupported difference scale", field="scale")
        percent_fraction = data.get("percent_fraction")
        if scale == "percent":
            if percent_fraction is None:
                raise InputError("percent scale requires percent_fraction")
            percent_fraction = number(
                percent_fraction, "percent_fraction", positive=True
            )
        replicates = integer(
            data.get("bootstrap_replicates", 2000),
            "bootstrap_replicates",
            MIN_BOOTSTRAP,
            MAX_BOOTSTRAP,
        )
        seed = integer(data.get("bootstrap_seed", 0), "bootstrap_seed", 0, 2**53 - 1)
        ratio = data.get("variance_ratio")
        if ratio is not None:
            ratio = number(ratio, "variance_ratio", positive=True)
        target = data.get("target_unit")
        if target is not None:
            identifier(target, "target_unit")
        return cls(
            identifier(data["id"]),
            reference[0],
            reference[1],
            candidate[0],
            candidate[1],
            start,
            end,
            specimen_tag,
            strategy,
            pairing_seconds,
            acceptable,
            range_low,
            range_high,
            scale,
            percent_fraction,
            replicates,
            seed,
            ratio,
            target,
        )

    def as_dict(self):
        return {
            "id": self.id,
            "reference": {
                "instrument": self.reference_instrument,
                "channel": self.reference_channel,
            },
            "candidate": {
                "instrument": self.candidate_instrument,
                "channel": self.candidate_channel,
            },
            "start": time_text(self.start),
            "end": time_text(self.end),
            "specimen_tag": self.specimen_tag,
            "replicate_strategy": self.replicate_strategy,
            "pairing_seconds": self.pairing_seconds,
            "acceptable_difference": decimal_text(self.acceptable_difference),
            "measurement_range": {
                "low": decimal_text(self.range_low),
                "high": decimal_text(self.range_high),
            },
            "scale": self.scale,
            "percent_fraction": (
                decimal_text(self.percent_fraction)
                if self.percent_fraction is not None
                else None
            ),
            "bootstrap_replicates": self.bootstrap_replicates,
            "bootstrap_seed": self.bootstrap_seed,
            "variance_ratio": (
                decimal_text(self.variance_ratio)
                if self.variance_ratio is not None
                else None
            ),
            "target_unit": self.target_unit,
        }


def _method(value, name):
    require_object(value, {"instrument", "channel"}, {"instrument", "channel"})
    return identifier(value["instrument"], name), identifier(value["channel"], name)


class PairingEngine:
    """Builds keyed, time-constrained pairs from accepted measurements."""

    def __init__(self, reports, catalog):
        self.reports = reports
        self.catalog = catalog

    def _load(self, instrument, channel, start, end):
        cursor = None
        items = []
        while True:
            page = self.reports.query(
                instrument=instrument,
                channel=channel,
                start=time_text(start),
                end=time_text(end),
                after=cursor,
                limit=10000,
            )
            for row in page["items"]:
                items.append(row)
            cursor = page["next"]
            if cursor is None:
                break
        if len(items) > MAX_POINTS:
            raise InputError("study supports at most ten thousand readings per method")
        return items

    def pair(self, spec):
        registry = self.catalog.registry()
        # Confirm both channels exist before reading anything.
        ref_unit = self.catalog.channel(
            spec.reference_instrument, spec.reference_channel
        )
        cand_unit = self.catalog.channel(
            spec.candidate_instrument, spec.candidate_channel
        )
        target_unit = spec.target_unit or ref_unit
        registry.get(target_unit)
        registry.validate_pair(ref_unit, target_unit)
        registry.validate_pair(cand_unit, target_unit)

        reference_rows = self._load(
            spec.reference_instrument, spec.reference_channel, spec.start, spec.end
        )
        candidate_rows = self._load(
            spec.candidate_instrument, spec.candidate_channel, spec.start, spec.end
        )
        exclusions = []
        ref_accepted, cand_accepted = self._accepted(
            reference_rows, "reference", exclusions, registry, target_unit
        ), self._accepted(
            candidate_rows, "candidate", exclusions, registry, target_unit
        )
        ref_eligible = self._check_certificates(
            ref_accepted, "reference", exclusions
        )
        cand_eligible = self._check_certificates(
            cand_accepted, "candidate", exclusions
        )

        ref_by_key = self._group_by_specimen(
            ref_eligible, spec, exclusions, "reference"
        )
        cand_by_key = self._group_by_specimen(
            cand_eligible, spec, exclusions, "candidate"
        )

        pairs = []
        all_keys = sorted(set(ref_by_key) | set(cand_by_key))
        for key in all_keys:
            ref_group = ref_by_key.get(key)
            cand_group = cand_by_key.get(key)
            if ref_group is None:
                for item in cand_group:
                    exclusions.append(
                        _exclusion(
                            "unpaired_specimen",
                            "candidate",
                            item,
                            {"specimen": key},
                        )
                    )
                continue
            if cand_group is None:
                for item in ref_group:
                    exclusions.append(
                        _exclusion(
                            "unpaired_specimen",
                            "reference",
                            item,
                            {"specimen": key},
                        )
                    )
                continue
            matched = self._match_within_specimen(ref_group, cand_group, spec, key)
            pairs.extend(matched["pairs"])
            exclusions.extend(matched["exclusions"])
        pairs.sort(key=lambda point: (point.specimen, point.replicate))
        return {
            "unit": target_unit,
            "points": pairs,
            "exclusions": exclusions,
            "reference_seen": len(reference_rows),
            "candidate_seen": len(candidate_rows),
            "reference_accepted": len(ref_accepted),
            "candidate_accepted": len(cand_accepted),
            "reference_eligible": len(ref_eligible),
            "candidate_eligible": len(cand_eligible),
        }

    def _accepted(self, rows, side, exclusions, registry, target_unit):
        accepted = []
        for row in rows:
            if row["measurement"] is None:
                exclusions.append(
                    {
                        "code": "processing_failure",
                        "side": side,
                        "observation": row["id"],
                        "message": "observation has no processed measurement",
                        "evidence": {"error": row["error"]},
                    }
                )
                continue
            measurement = Measurement.restore(row["measurement"])
            if not measurement.accepted:
                exclusions.append(
                    {
                        "code": "rejected_measurement",
                        "side": side,
                        "observation": row["id"],
                        "message": "measurement failed the quality policy",
                        "evidence": {
                            "issues": [issue.as_dict() for issue in measurement.issues]
                        },
                    }
                )
                continue
            converted = registry.convert(
                measurement.value, measurement.unit, target_unit
            )
            accepted.append((row, measurement, converted))
        return accepted

    def _check_certificates(self, accepted, side, exclusions):
        book = self.catalog.certificates()
        active = {cert.id for cert in book.active()}
        eligible = []
        for triple in accepted:
            _row, measurement, _converted = triple
            if measurement.certificate not in active:
                exclusions.append(
                    {
                        "code": "missing_or_revoked_certificate",
                        "side": side,
                        "observation": triple[0]["id"],
                        "message": "calibration certificate is absent or revoked",
                        "evidence": {"certificate": measurement.certificate},
                    }
                )
                continue
            eligible.append(triple)
        return eligible

    def _group_by_specimen(self, accepted, spec, exclusions, side):
        grouped = {}
        for row, measurement, converted in accepted:
            key = measurement.observation.tags.get(spec.specimen_tag)
            if key is None:
                exclusions.append(
                    {
                        "code": "missing_specimen_key",
                        "side": side,
                        "observation": row["id"],
                        "message": "observation tag for the specimen key is absent",
                        "evidence": {"specimen_tag": spec.specimen_tag},
                    }
                )
                continue
            grouped.setdefault(key, []).append((row, measurement, converted))
        return grouped

    def _match_within_specimen(self, ref_group, cand_group, spec, specimen):
        # Greedy globally-closest matching; equidistant ties are reported,
        # never guessed, and the readings involved stay unmatched.
        edges = []
        for ri, (r_row, _r_m, _r_v) in enumerate(ref_group):
            r_time = instant(r_row["event_time"])
            for ci, (c_row, _c_m, _c_v) in enumerate(cand_group):
                seconds = abs(
                    (instant(c_row["event_time"]) - r_time).total_seconds()
                )
                if seconds <= spec.pairing_seconds:
                    edges.append((seconds, ri, ci))
        edges.sort(key=lambda edge: (edge[0], edge[1], edge[2]))
        ambiguous_r = set()
        exclusions = []
        for ri in range(len(ref_group)):
            candidates = [edge for edge in edges if edge[1] == ri]
            if not candidates:
                continue
            closest = candidates[0][0]
            tied = [edge for edge in candidates if edge[0] == closest]
            if len(tied) > 1:
                ambiguous_r.add(ri)
                exclusions.append(
                    {
                        "code": "ambiguous_pair",
                        "side": "reference",
                        "observation": ref_group[ri][0]["id"],
                        "message": "multiple candidate readings at equal time distance",
                        "evidence": {
                            "specimen": specimen,
                            "candidates": [
                                cand_group[ci][0]["id"]
                                for _seconds, _r, ci in tied
                            ],
                        },
                    }
                )
        used_r, used_c = set(), set()
        pairs = []
        for seconds, ri, ci in edges:
            if ri in used_r or ci in used_c or ri in ambiguous_r:
                continue
            r_row, r_measurement, r_value = ref_group[ri]
            c_row, c_measurement, c_value = cand_group[ci]
            used_r.add(ri)
            used_c.add(ci)
            pairs.append(
                Point(
                    specimen=specimen,
                    replicate=str(len(pairs) + 1),
                    reference=r_value,
                    candidate=c_value,
                    reference_id=r_row["id"],
                    candidate_id=c_row["id"],
                    reference_revision=r_row["revision"],
                    candidate_revision=c_row["revision"],
                    certificate_reference=r_measurement.certificate,
                    certificate_candidate=c_measurement.certificate,
                    reference_at=r_row["event_time"],
                    candidate_at=c_row["event_time"],
                )
            )
        for ri, triple in enumerate(ref_group):
            if ri not in used_r and ri not in ambiguous_r:
                exclusions.append(
                    _exclusion(
                        "pair_outside_time_window",
                        "reference",
                        triple,
                        {"specimen": specimen, "pairing_seconds": spec.pairing_seconds},
                    )
                )
        for ci, triple in enumerate(cand_group):
            if ci not in used_c:
                exclusions.append(
                    _exclusion(
                        "pair_outside_time_window",
                        "candidate",
                        triple,
                        {"specimen": specimen, "pairing_seconds": spec.pairing_seconds},
                    )
                )
        return {"pairs": pairs, "exclusions": exclusions}


def _exclusion(code, side, triple, evidence):
    row, _measurement, _value = triple
    return {
        "code": code,
        "side": side,
        "observation": row["id"],
        "message": code,
        "evidence": evidence,
    }


class AgreementStudies:
    def __init__(self, database, reports, catalog):
        self.db = database
        self.reports = reports
        self.catalog = catalog
        with self.db.transaction() as con:
            con.execute(
                """CREATE TABLE IF NOT EXISTS agreement_studies(id TEXT PRIMARY KEY,
                specification TEXT NOT NULL, fingerprint TEXT NOT NULL, result TEXT NOT NULL,
                version INTEGER NOT NULL)"""
            )

    def run(self, data):
        spec = StudySpec.parse(data)
        pairing = PairingEngine(self.reports, self.catalog).pair(spec)
        result = self._analyze(spec, pairing)
        body = {
            "specification": spec.as_dict(),
            "pairing": {
                "unit": pairing["unit"],
                "reference_seen": pairing["reference_seen"],
                "candidate_seen": pairing["candidate_seen"],
                "reference_accepted": pairing["reference_accepted"],
                "candidate_accepted": pairing["candidate_accepted"],
                "reference_eligible": pairing["reference_eligible"],
                "candidate_eligible": pairing["candidate_eligible"],
                "paired_readings": len(pairing["points"]),
                "paired_specimens": len({p.specimen for p in pairing["points"]}),
                "exclusions": pairing["exclusions"],
            },
            "analysis": result,
        }
        fingerprint = digest(body)
        with self.db.transaction() as con:
            previous = con.execute(
                "SELECT fingerprint,version FROM agreement_studies WHERE id=?",
                (spec.id,),
            ).fetchone()
            if previous:
                if previous["fingerprint"] != fingerprint:
                    raise Conflict("agreement study id is immutable; create a new version id")
                return self.get(spec.id)
            con.execute(
                "INSERT INTO agreement_studies VALUES(?,?,?,?,1)",
                (
                    spec.id,
                    canonical(spec.as_dict()),
                    fingerprint,
                    canonical(body),
                ),
            )
            append(
                con,
                "agreement_study_run",
                spec.id,
                {
                    "fingerprint": fingerprint,
                    "specimens": body["pairing"]["paired_specimens"],
                    "decision": result["decision"]["status"],
                },
            )
        return self.get(spec.id)

    def _analyze(self, spec, pairing):
        points = pairing["points"]
        diagnostics = list(pairing_to_diagnostics(pairing))
        if not points:
            return {
                "scale": spec.scale,
                "bland_altman": None,
                "deming": None,
                "bootstrap": None,
                "decision": {
                    "status": "inconclusive",
                    "reasons": ["no_paired_readings"],
                },
                "diagnostics": diagnostics,
            }
        try:
            transformed = transform(
                points, spec.scale, spec.percent_fraction
            )
        except InputError as exc:
            diagnostics.append(exc.as_dict())
            return {
                "scale": spec.scale,
                "bland_altman": None,
                "deming": None,
                "bootstrap": None,
                "decision": {
                    "status": "inconclusive",
                    "reasons": ["scale_not_applicable"],
                },
                "diagnostics": diagnostics,
            }
        diagnostics.extend(transformed.diagnostics)
        summary = bland_altman(transformed)
        # Specimen-level aggregation for independent-unit statistics.
        specimen_points = aggregate_specimens(points, spec.replicate_strategy)
        specimen_transformed = transform(
            specimen_points, spec.scale, spec.percent_fraction
        )
        specimen_summary = bland_altman(specimen_transformed)
        bootstrap = None
        reasons = []
        if len(specimen_points) >= 2:
            try:
                bootstrap = cluster_bootstrap(
                    specimen_transformed,
                    ba_bootstrap_statistic(spec.scale, spec.percent_fraction),
                    replicates=spec.bootstrap_replicates,
                    seed=spec.bootstrap_seed,
                )
            except InputError as exc:
                diagnostics.append(exc.as_dict())
        else:
            reasons.append("insufficient_specimens_for_bootstrap")

        deming_result = None
        if spec.variance_ratio is not None:
            if len({p.reference for p in specimen_points}) < 2:
                reasons.append("degenerate_deming_data")
                diagnostics.append(
                    {
                        "code": "degenerate_deming_data",
                        "message": "Deming regression needs at least two distinct values",
                    }
                )
            else:
                try:
                    deming_result = deming(
                        [p.reference for p in specimen_points],
                        [p.candidate for p in specimen_points],
                        spec.variance_ratio,
                    )
                except InputError as exc:
                    reasons.append("deming_failed")
                    diagnostics.append(exc.as_dict())

        coverage = range_coverage(specimen_points, spec.range_low, spec.range_high)
        decision = decide(
            spec,
            specimen_summary,
            bootstrap,
            coverage,
            len(specimen_points),
            reasons,
        )
        return {
            "scale": spec.scale,
            "reading_level": summary,
            "bland_altman": specimen_summary,
            "deming": deming_result,
            "bootstrap": bootstrap,
            "range_coverage": coverage,
            "decision": decision,
            "diagnostics": diagnostics,
            "points": [point_as_dict(point) for point in points],
        }

    def get(self, study_id):
        identifier(study_id)
        row = self.db.one(
            "SELECT * FROM agreement_studies WHERE id=?", (study_id,)
        )
        if row is None:
            raise Missing("agreement study not found")
        payload = strict_json(row["result"])
        if digest(payload) != row["fingerprint"]:
            raise StateError("stored agreement study fingerprint mismatch")
        return {
            "id": study_id,
            "version": row["version"],
            "fingerprint": row["fingerprint"],
            "payload": payload,
        }

    def list(self):
        return [
            {
                "id": row["id"],
                "version": row["version"],
                "fingerprint": row["fingerprint"],
            }
            for row in self.db.rows(
                "SELECT id,version,fingerprint FROM agreement_studies ORDER BY id"
            )
        ]

    def compare_versions(self, left_id, right_id):
        left = self.get(left_id)["payload"]
        right = self.get(right_id)["payload"]
        changed_specimens = []
        left_points = {p["specimen"]: p for p in left["analysis"].get("points", [])}
        right_points = {p["specimen"]: p for p in right["analysis"].get("points", [])}
        for key in sorted(set(left_points) | set(right_points)):
            if left_points.get(key) != right_points.get(key):
                changed_specimens.append(key)
        return {
            "left": left_id,
            "right": right_id,
            "left_decision": left["analysis"]["decision"],
            "right_decision": right["analysis"]["decision"],
            "decision_changed": (
                left["analysis"]["decision"]["status"]
                != right["analysis"]["decision"]["status"]
            ),
            "specimens_added": sorted(set(right_points) - set(left_points)),
            "specimens_removed": sorted(set(left_points) - set(right_points)),
            "specimens_changed": changed_specimens,
        }


def aggregate_specimens(points, strategy):
    groups = {}
    for point in points:
        groups.setdefault(point.specimen, []).append(point)
    result = []
    for specimen, group in sorted(groups.items()):
        ordered = sorted(group, key=lambda p: (p.reference_at, p.reference_id))
        if strategy == "first":
            chosen = ordered[0]
        elif strategy == "last":
            chosen = ordered[-1]
        else:
            with localcontext() as ctx:
                ctx.prec = 40
                ref = sum((p.reference for p in ordered), number("0")) / number(
                    str(len(ordered))
                )
                cand = sum((p.candidate for p in ordered), number("0")) / number(
                    str(len(ordered))
                )
            chosen = Point(
                specimen=specimen,
                replicate="mean",
                reference=ref,
                candidate=cand,
                reference_id=";".join(p.reference_id for p in ordered),
                candidate_id=";".join(p.candidate_id for p in ordered),
                reference_revision=ordered[0].reference_revision,
                candidate_revision=ordered[0].candidate_revision,
                certificate_reference=ordered[0].certificate_reference,
                certificate_candidate=ordered[0].certificate_candidate,
                reference_at=ordered[0].reference_at,
                candidate_at=ordered[0].candidate_at,
            )
        result.append(
            Point(
                specimen=specimen,
                replicate=chosen.replicate,
                reference=chosen.reference,
                candidate=chosen.candidate,
                reference_id=chosen.reference_id,
                candidate_id=chosen.candidate_id,
                reference_revision=chosen.reference_revision,
                candidate_revision=chosen.candidate_revision,
                certificate_reference=chosen.certificate_reference,
                certificate_candidate=chosen.certificate_candidate,
                reference_at=chosen.reference_at,
                candidate_at=chosen.candidate_at,
            )
        )
    return result


def range_coverage(points, low, high):
    if not points:
        return {
            "low": decimal_text(low),
            "high": decimal_text(high),
            "observed_low": None,
            "observed_high": None,
            "specimens_inside": 0,
            "coverage_fraction": "0",
            "fully_covered": False,
        }
    observed_low = min(p.reference for p in points)
    observed_high = max(p.reference for p in points)
    inside = [p for p in points if low <= p.reference <= high]
    return {
        "low": decimal_text(low),
        "high": decimal_text(high),
        "observed_low": decimal_text(observed_low),
        "observed_high": decimal_text(observed_high),
        "specimens_inside": len(inside),
        "coverage_fraction": decimal_text(
            number(str(len(inside))) / number(str(len(points)))
        ),
        "fully_covered": observed_low <= low and observed_high >= high,
    }


def decide(spec, summary, bootstrap, coverage, specimen_count, reasons):
    reasons = list(reasons)
    if summary is None:
        return {
            "status": "inconclusive",
            "reasons": reasons or ["no_paired_readings"],
            "acceptable_difference": decimal_text(spec.acceptable_difference),
        }
    if specimen_count < 8:
        reasons.append("small_sample")
    if not coverage["fully_covered"]:
        reasons.append("range_not_covered")
    lower = number(summary["lower_limit"])
    upper = number(summary["upper_limit"])
    acceptable = spec.acceptable_difference
    pointwise_ok = -acceptable <= lower and upper <= acceptable
    interval_ok = pointwise_ok
    if bootstrap is not None and "lower_limit" in bootstrap and "upper_limit" in bootstrap:
        interval_ok = (
            number(bootstrap["lower_limit"]["lower"]) >= -acceptable
            and number(bootstrap["upper_limit"]["upper"]) <= acceptable
        )
    else:
        reasons.append("bootstrap_unavailable")
    if pointwise_ok and interval_ok and not reasons:
        status = "interchangeable"
    elif pointwise_ok and interval_ok:
        # Limits satisfy the allowance, but the evidence base is not yet
        # sufficient for an unconditional replacement claim.
        status = "conditionally_interchangeable"
    else:
        # A limit (or its interval) crosses the allowance: no qualification
        # level can support replacement regardless of other reasons.
        status = "not_interchangeable"
        reasons = [reason for reason in reasons if reason != "small_sample"]
    return {
        "status": status,
        "reasons": reasons,
        "acceptable_difference": decimal_text(acceptable),
        "limits_within_allowance": pointwise_ok,
        "bootstrap_limits_within_allowance": interval_ok,
        "lower_limit": summary["lower_limit"],
        "upper_limit": summary["upper_limit"],
        "bootstrap_intervals": bootstrap,
    }


def pairing_to_diagnostics(pairing):
    counts = {}
    for exclusion in pairing["exclusions"]:
        counts[exclusion["code"]] = counts.get(exclusion["code"], 0) + 1
    return (
        {
            "code": "exclusion_summary",
            "message": "pairing exclusions by reason",
            "counts": counts,
        },
    )


def point_as_dict(point):
    return {
        "specimen": point.specimen,
        "replicate": point.replicate,
        "reference": decimal_text(point.reference),
        "candidate": decimal_text(point.candidate),
        "difference": decimal_text(point.difference),
        "average": decimal_text(point.average),
        "reference_observation": point.reference_id,
        "candidate_observation": point.candidate_id,
        "reference_revision": point.reference_revision,
        "candidate_revision": point.candidate_revision,
        "reference_certificate": point.certificate_reference,
        "candidate_certificate": point.certificate_candidate,
        "reference_at": point.reference_at,
        "candidate_at": point.candidate_at,
    }
