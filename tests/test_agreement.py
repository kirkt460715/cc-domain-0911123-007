import json
import math
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

from metrology.agreement.statistics import (
    Point,
    bland_altman,
    cluster_bootstrap,
    deming,
    transform,
)
from metrology.agreement.studies import (
    AgreementStudies,
    StudySpec,
    aggregate_specimens,
)
from metrology.agreement.export import AgreementExports
from metrology.console import dispatch
from metrology.errors import Conflict, InputError, Missing
from support import (
    RECEIVED,
    START,
    application,
    certificate,
    observation,
    rows,
)

END = "2026-11-11T00:02:00Z"


def portable_setup(app, cert_id="pcert"):
    app.catalog.add_instrument(
        {"id": "portable", "label": "Portable", "channels": {"p": "kPa"}}
    )
    app.catalog.add_source("psrc")
    app.catalog.add_certificate(
        certificate(id=cert_id, instrument="portable")
    )
    app.catalog.bind("portable", "p", "q1")


def pair(seq, ref_id, port_id, base, port_value, key, minute_offset=0):
    t = f"2026-11-11T00:00:{minute_offset:02d}Z"
    return [
        observation(
            seq,
            id=ref_id,
            value=str(base),
            event_time=t,
            tags={"lot": "A", "specimen": key},
        ),
        observation(
            seq + 1,
            id=port_id,
            instrument="portable",
            source="psrc",
            value=str(port_value),
            event_time=t,
            tags={"lot": "A", "specimen": key},
        ),
    ]


def base_spec(study_id, **changes):
    spec = {
        "id": study_id,
        "reference": {"instrument": "i1", "channel": "p"},
        "candidate": {"instrument": "portable", "channel": "p"},
        "start": START,
        "end": END,
        "specimen_tag": "specimen",
        "replicate_strategy": "mean",
        "pairing_seconds": 2,
        "acceptable_difference": "3",
        "measurement_range": {"low": "5", "high": "16"},
        "scale": "absolute",
        "bootstrap_replicates": 100,
        "bootstrap_seed": 7,
    }
    spec.update(changes)
    return spec


def ingest_grid(app, offset, count, port_function):
    data = []
    seq = offset
    for i in range(count):
        base = Decimal(5 + i)
        data.extend(
            pair(
                seq,
                f"b{i}",
                f"p{i}",
                base,
                port_function(base, i),
                f"S{i:02d}",
                minute_offset=i,
            )
        )
        seq += 2
    app.pipeline.ingest("b", rows(*data), received_at=RECEIVED)


class FixedBiasTests(unittest.TestCase):
    """High correlation with a constant offset must not pass on correlation."""

    def setUp(self):
        self.app = application()
        self.addCleanup(self.app.close)
        portable_setup(self.app)
        ingest_grid(self.app, 0, 12, lambda base, i: base + 2)
        self.studies = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        )

    def test_bland_altman_reveals_bias_despite_perfect_correlation(self):
        result = self.studies.run(base_spec("fixed"))["payload"]
        ba = result["analysis"]["bland_altman"]
        self.assertEqual(ba["correlation"], "1")
        self.assertEqual(ba["bias"], "2")
        self.assertEqual(ba["lower_limit"], "2")
        self.assertEqual(ba["upper_limit"], "2")
        self.assertEqual(ba["proportional_bias"]["slope"], "0")
        self.assertEqual(ba["proportional_bias"]["intercept"], "2")
        self.assertEqual(result["pairing"]["paired_specimens"], 12)
        self.assertEqual(result["pairing"]["exclusions"], [])

    def test_decision_compares_limits_with_allowance(self):
        inside = self.studies.run(base_spec("ok"))["payload"]
        self.assertEqual(
            inside["analysis"]["decision"]["status"], "interchangeable"
        )
        outside = self.studies.run(
            base_spec("tight", acceptable_difference="1.5")
        )["payload"]
        self.assertEqual(
            outside["analysis"]["decision"]["status"], "not_interchangeable"
        )
        self.assertFalse(
            outside["analysis"]["decision"]["limits_within_allowance"]
        )

    def test_deming_recovers_offset_without_ordinary_least_squares(self):
        result = self.studies.run(
            base_spec("deming", variance_ratio="1")
        )["payload"]
        deming_result = result["analysis"]["deming"]
        self.assertEqual(deming_result["slope"], "1")
        self.assertEqual(deming_result["intercept"], "2")
        self.assertEqual(deming_result["variance_ratio"], "1")
        self.assertIn("slope_interval", deming_result)

    def test_bootstrap_is_specimen_level_and_deterministic(self):
        first = self.studies.run(base_spec("boot"))["payload"]
        second = self.studies.run(base_spec("boot2"))["payload"]
        self.assertEqual(
            first["analysis"]["bootstrap"], second["analysis"]["bootstrap"]
        )
        intervals = first["analysis"]["bootstrap"]
        self.assertEqual(intervals["bias"]["lower"], "2")
        self.assertEqual(intervals["bias"]["replicates"], 100)


class ProportionalBiasTests(unittest.TestCase):
    def setUp(self):
        self.app = application()
        self.addCleanup(self.app.close)
        portable_setup(self.app)
        ingest_grid(
            self.app,
            0,
            10,
            lambda base, i: (base * Decimal("1.2")).quantize(Decimal("0.001")),
        )
        self.studies = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        )

    def test_proportional_bias_fails_even_with_high_correlation(self):
        result = self.studies.run(
            base_spec(
                "prop",
                acceptable_difference="1",
                measurement_range={"low": "5", "high": "14"},
                variance_ratio="2",
            )
        )["payload"]
        ba = result["analysis"]["bland_altman"]
        self.assertNotEqual(ba["proportional_bias"]["slope"], "0")
        slope = Decimal(result["analysis"]["deming"]["slope"])
        self.assertLess(slope, Decimal("1"))
        self.assertEqual(
            result["analysis"]["decision"]["status"], "not_interchangeable"
        )

    def test_percent_and_log_scales(self):
        percent = self.studies.run(
            base_spec(
                "pct",
                scale="percent",
                percent_fraction="100",
                acceptable_difference="50",
            )
        )["payload"]
        self.assertEqual(
            percent["analysis"]["bland_altman"]["scale"], "percent"
        )
        log = self.studies.run(
            base_spec("log", scale="log", acceptable_difference="1")
        )["payload"]
        self.assertEqual(log["analysis"]["bland_altman"]["scale"], "log")
        diagnostic_codes = {
            item["code"] for item in log["analysis"]["diagnostics"]
        }
        self.assertIn("log_scale", diagnostic_codes)


class RepeatedMeasurementTests(unittest.TestCase):
    def setUp(self):
        self.app = application()
        self.addCleanup(self.app.close)
        portable_setup(self.app)
        data = []
        seq = 0
        for i in range(10):
            base = Decimal(5 + i)
            for replicate in range(2):
                t = f"2026-11-11T00:00:{i * 2 + replicate:02d}Z"
                noise = Decimal("0.1") if replicate else Decimal("-0.1")
                data.extend(
                    pair(
                        seq,
                        f"b{i}-{replicate}",
                        f"p{i}-{replicate}",
                        base,
                        (base * Decimal("1.2") + noise).quantize(
                            Decimal("0.001")
                        ),
                        f"S{i:02d}",
                        minute_offset=i * 2 + replicate,
                    )
                )
                seq += 2
        self.app.pipeline.ingest("b", rows(*data), received_at=RECEIVED)
        self.studies = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        )

    def test_repeats_are_one_independent_specimen(self):
        result = self.studies.run(
            base_spec(
                "rep",
                acceptable_difference="1",
                measurement_range={"low": "5", "high": "14"},
            )
        )["payload"]
        self.assertEqual(result["pairing"]["paired_readings"], 20)
        self.assertEqual(result["pairing"]["paired_specimens"], 10)
        reading = result["analysis"]["reading_level"]
        self.assertEqual(reading["specimens"], 10)
        self.assertEqual(reading["repeated_specimens"], 10)
        self.assertEqual(reading["repeatability_df"], 10)
        self.assertIsNotNone(reading["repeatability_limit"])
        specimen = result["analysis"]["bland_altman"]
        self.assertEqual(specimen["pairs"], 10)
        # mean strategy averages the ±0.1 noise away
        self.assertEqual(
            specimen["components"]["within_specimen_variance"], "0"
        )

    def test_replicate_strategies_choose_distinct_readings(self):
        from metrology.agreement.studies import PairingEngine

        spec = StudySpec.parse(base_spec("raw"))
        points = PairingEngine(self.app.reports, self.app.catalog).pair(spec)[
            "points"
        ]
        first = aggregate_specimens(points, "first")[0]
        last = aggregate_specimens(points, "last")[0]
        mean = aggregate_specimens(points, "mean")[0]
        self.assertIn("b0-0", first.reference_id)
        self.assertIn("b0-1", last.reference_id)
        self.assertEqual(mean.reference_id, "b0-0;b0-1")
        self.assertEqual(mean.replicate, "mean")


class PairingDiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.app = application()
        self.addCleanup(self.app.close)
        portable_setup(self.app)
        data = []
        seq = 0
        for i in range(9):
            data.extend(
                pair(
                    seq,
                    f"b{i}",
                    f"p{i}",
                    Decimal(5 + i),
                    Decimal(5.5 + i),
                    f"S{i:02d}",
                    minute_offset=i,
                )
            )
            seq += 2
        # reference-only specimen
        data.append(
            observation(
                seq,
                id="b99",
                value="8",
                event_time="2026-11-11T00:00:30Z",
                tags={"lot": "A", "specimen": "S99"},
            )
        )
        # candidate outside the time window for specimen S00
        data.append(
            observation(
                seq + 1,
                id="p00-late",
                instrument="portable",
                source="psrc",
                value="6.5",
                event_time="2026-11-11T00:00:55Z",
                tags={"lot": "A", "specimen": "S00"},
            )
        )
        # missing specimen tag
        data.append(
            observation(
                seq + 2,
                id="pnokey",
                instrument="portable",
                source="psrc",
                value="6.5",
                event_time="2026-11-11T00:00:05Z",
                tags={"lot": "A"},
            )
        )
        self.app.pipeline.ingest("b", rows(*data), received_at=RECEIVED)
        self.studies = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        )

    def test_exclusions_are_itemized_without_positional_matching(self):
        result = self.studies.run(base_spec("diag"))["payload"]
        codes = {}
        for exclusion in result["pairing"]["exclusions"]:
            codes[exclusion["code"]] = codes.get(exclusion["code"], 0) + 1
        self.assertEqual(codes["unpaired_specimen"], 1)
        self.assertEqual(codes["pair_outside_time_window"], 1)
        self.assertEqual(codes["missing_specimen_key"], 1)
        self.assertEqual(result["pairing"]["paired_specimens"], 9)

    def test_single_observation_is_never_consumed_twice(self):
        app = application()
        self.addCleanup(app.close)
        portable_setup(app)
        app.pipeline.ingest(
            "b",
            rows(
                observation(
                    0,
                    id="b1",
                    value="5",
                    event_time="2026-11-11T00:00:00Z",
                    tags={"lot": "A", "specimen": "S1"},
                ),
                observation(
                    1,
                    id="b2",
                    value="6",
                    event_time="2026-11-11T00:00:01Z",
                    tags={"lot": "A", "specimen": "S1"},
                ),
                observation(
                    2,
                    id="p1",
                    instrument="portable",
                    source="psrc",
                    value="5.5",
                    event_time="2026-11-11T00:00:01Z",
                    tags={"lot": "A", "specimen": "S1"},
                ),
            ),
            received_at=RECEIVED,
        )
        result = AgreementStudies(
            app.database, app.reports, app.catalog
        ).run(base_spec("once", pairing_seconds=5))["payload"]
        self.assertEqual(result["pairing"]["paired_readings"], 1)
        self.assertEqual(
            [item["code"] for item in result["pairing"]["exclusions"]],
            ["pair_outside_time_window"],
        )

    def test_equidistant_candidates_are_ambiguous_not_guessed(self):
        app = application()
        self.addCleanup(app.close)
        portable_setup(app)
        app.pipeline.ingest(
            "b",
            rows(
                observation(
                    0,
                    id="b1",
                    value="5",
                    event_time="2026-11-11T00:00:00Z",
                    tags={"lot": "A", "specimen": "S1"},
                ),
                observation(
                    1,
                    id="p1",
                    instrument="portable",
                    source="psrc",
                    value="5.1",
                    event_time="2026-11-11T00:00:01Z",
                    tags={"lot": "A", "specimen": "S1"},
                ),
                observation(
                    2,
                    id="p2",
                    instrument="portable",
                    source="psrc",
                    value="5.2",
                    event_time="2026-11-11T00:00:01Z",
                    tags={"lot": "A", "specimen": "S1"},
                ),
            ),
            received_at=RECEIVED,
        )
        result = AgreementStudies(
            app.database, app.reports, app.catalog
        ).run(base_spec("amb", pairing_seconds=5))["payload"]
        self.assertEqual(result["pairing"]["paired_readings"], 0)
        self.assertEqual(
            [item["code"] for item in result["pairing"]["exclusions"]][0],
            "ambiguous_pair",
        )


class RevokedCertificateTests(unittest.TestCase):
    def test_revoked_certificate_readings_are_excluded(self):
        app = application()
        self.addCleanup(app.close)
        portable_setup(app)
        ingest_grid(app, 0, 9, lambda base, i: base + Decimal("0.5"))
        app.catalog.revoke_certificate("pcert", "recall")
        result = AgreementStudies(
            app.database, app.reports, app.catalog
        ).run(base_spec("rev"))["payload"]
        self.assertEqual(result["pairing"]["paired_readings"], 0)
        codes = {item["code"] for item in result["pairing"]["exclusions"]}
        self.assertIn("missing_or_revoked_certificate", codes)
        self.assertEqual(
            result["analysis"]["decision"]["status"], "inconclusive"
        )


class SmallSampleAndCoverageTests(unittest.TestCase):
    def setUp(self):
        self.app = application()
        self.addCleanup(self.app.close)
        portable_setup(self.app)
        ingest_grid(self.app, 0, 5, lambda base, i: base + Decimal("0.2"))

    def test_small_sample_is_conditional_never_auto_pass(self):
        result = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        ).run(
            base_spec(
                "small",
                measurement_range={"low": "5", "high": "9"},
            )
        )["payload"]
        self.assertEqual(
            result["analysis"]["decision"]["status"],
            "conditionally_interchangeable",
        )
        self.assertIn("small_sample", result["analysis"]["decision"]["reasons"])

    def test_uncovered_range_blocks_interchangeability(self):
        result = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        ).run(
            base_spec(
                "coverage",
                measurement_range={"low": "100", "high": "200"},
            )
        )["payload"]
        self.assertFalse(result["analysis"]["range_coverage"]["fully_covered"])
        self.assertIn(
            "range_not_covered", result["analysis"]["decision"]["reasons"]
        )

    def test_small_sample_with_failing_limits_is_not_conditional_pass(self):
        # allowance far below the 0.2 bias: limits cross the allowance, so
        # the small sample must never be promoted to a conditional pass
        result = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        ).run(
            base_spec(
                "failing",
                acceptable_difference="0.05",
                measurement_range={"low": "5", "high": "9"},
            )
        )["payload"]
        self.assertEqual(
            result["analysis"]["decision"]["status"], "not_interchangeable"
        )
        self.assertNotIn(
            "small_sample", result["analysis"]["decision"]["reasons"]
        )


class ScaleDiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.app = application()
        self.addCleanup(self.app.close)
        portable_setup(self.app)
        data = []
        values = [(0, 1, "S00"), (-2, -1, "S01"), (5, 6, "S02"), (7, 8, "S03")]
        for index, (ref, cand, key) in enumerate(values):
            t = f"2026-11-11T00:00:{index:02d}Z"
            data.append(
                observation(
                    index,
                    id=f"b{index}",
                    value=str(ref),
                    event_time=t,
                    tags={"lot": "A", "specimen": key},
                )
            )
            data.append(
                observation(
                    index + 10,
                    id=f"p{index}",
                    instrument="portable",
                    source="psrc",
                    value=str(cand),
                    event_time=t,
                    tags={"lot": "A", "specimen": key},
                )
            )
        self.app.pipeline.ingest("b", rows(*data), received_at=RECEIVED)

    def test_percent_scale_reports_inconclusive_on_nonpositive_reference(self):
        result = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        ).run(
            base_spec("pct", scale="percent", percent_fraction="100")
        )["payload"]
        self.assertEqual(
            result["analysis"]["decision"]["status"], "inconclusive"
        )
        self.assertIn(
            "scale_not_applicable",
            result["analysis"]["decision"]["reasons"],
        )
        self.assertIsNone(result["analysis"]["bland_altman"])

    def test_log_scale_reports_inconclusive_with_diagnostic(self):
        result = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        ).run(base_spec("log", scale="log"))["payload"]
        self.assertEqual(
            result["analysis"]["decision"]["status"], "inconclusive"
        )
        self.assertIn(
            "scale_not_applicable",
            result["analysis"]["decision"]["reasons"],
        )
        self.assertIsNone(result["analysis"]["bland_altman"])


class DemingValidationTests(unittest.TestCase):
    def test_variance_ratio_must_be_positive(self):
        for ratio in ("0", "-1", "nonsense"):
            with self.assertRaises(InputError):
                StudySpec.parse(base_spec("bad", variance_ratio=ratio))

    def test_deming_rejects_degenerate_data(self):
        points = [
            Point(
                f"S{i}",
                "1",
                Decimal("5"),
                Decimal("5.1"),
                f"b{i}",
                f"p{i}",
                1,
                1,
                "c1",
                "pcert",
                START,
                START,
            )
            for i in range(9)
        ]
        with self.assertRaises(InputError):
            deming(
                [p.reference for p in points],
                [p.candidate for p in points],
                "1",
            )

    def test_deming_matches_closed_form_against_manual_calculation(self):
        # y = 2x + 1 with symmetric noise; lambda=1 closed form
        refs = [Decimal(v) for v in (1, 2, 3, 4, 5, 6)]
        cands = [Decimal(v) * 2 + 1 for v in refs]
        result = deming(refs, cands, "1")
        self.assertEqual(result["slope"], "2")
        self.assertEqual(result["intercept"], "1")


class StatisticsUnitTests(unittest.TestCase):
    def _points(self):
        return [
            Point(
                "S1",
                "1",
                Decimal("10"),
                Decimal("12"),
                "b1",
                "p1",
                1,
                1,
                "c1",
                "pc",
                START,
                START,
            ),
            Point(
                "S2",
                "1",
                Decimal("20"),
                Decimal("22"),
                "b2",
                "p2",
                1,
                1,
                "c1",
                "pc",
                START,
                START,
            ),
        ]

    def test_cluster_bootstrap_resamples_whole_specimens(self):
        transformed = transform(self._points(), "absolute")
        intervals = cluster_bootstrap(
            transformed,
            lambda sampled: {"bias": Decimal(
                sum(p.difference for p in sampled) / Decimal(len(sampled))
            )},
            replicates=50,
            seed=11,
        )
        # every resampled difference is 2, so the interval is degenerate
        self.assertEqual(intervals["bias"]["lower"], "2")
        self.assertEqual(intervals["bias"]["upper"], "2")

    def test_within_specimen_variance_uses_repeated_differences(self):
        points = self._points() + [
            Point(
                "S1",
                "2",
                Decimal("10"),
                Decimal("13"),
                "b1b",
                "p1b",
                1,
                1,
                "c1",
                "pc",
                START,
                START,
            )
        ]
        summary = bland_altman(transform(points, "absolute"))
        # within S1 differences are 2 and 3; with one degree of freedom the
        # pooled within-specimen variance is ((2-2.5)^2+(3-2.5)^2)/1 = 0.5
        self.assertEqual(
            summary["components"]["within_specimen_variance"], "0.5"
        )
        self.assertEqual(summary["repeatability_df"], 1)
        self.assertIsNotNone(summary["repeatability_limit"])

    def test_log_transform_is_natural_log(self):
        transformed = transform(self._points(), "log")
        self.assertAlmostEqual(
            float(transformed.differences[0]),
            math.log(12 / 10),
            places=12,
        )


class ImmutabilityAndVersionTests(unittest.TestCase):
    def setUp(self):
        self.app = application()
        self.addCleanup(self.app.close)
        portable_setup(self.app)
        ingest_grid(self.app, 0, 9, lambda base, i: base + Decimal("0.3"))

    def test_study_id_is_immutable(self):
        studies = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        )
        first = studies.run(base_spec("v1"))
        reused = studies.run(base_spec("v1"))
        self.assertEqual(reused["fingerprint"], first["fingerprint"])
        with self.assertRaises(Conflict):
            studies.run(base_spec("v1", acceptable_difference="9"))

    def test_version_compare_lists_changed_specimens(self):
        studies = AgreementStudies(
            self.app.database, self.app.reports, self.app.catalog
        )
        studies.run(base_spec("v1", pairing_seconds=2))
        studies.run(base_spec("v2", pairing_seconds=60))
        comparison = studies.compare_versions("v1", "v2")
        self.assertFalse(comparison["decision_changed"])
        self.assertEqual(comparison["specimens_changed"], [])
        with self.assertRaises(Missing):
            studies.compare_versions("v1", "missing")

    def test_specification_validation(self):
        with self.assertRaises(InputError):
            StudySpec.parse(base_spec("bad", replicate_strategy="median"))
        with self.assertRaises(InputError):
            StudySpec.parse(base_spec("bad", scale="rank"))
        with self.assertRaises(InputError):
            StudySpec.parse(
                base_spec("bad", measurement_range={"low": "9", "high": "1"})
            )
        with self.assertRaises(InputError):
            StudySpec.parse(
                base_spec(
                    "same",
                    candidate={"instrument": "i1", "channel": "p"},
                )
            )


class ExportAndConsoleTests(unittest.TestCase):
    def setUp(self):
        self.app = application()
        self.addCleanup(self.app.close)
        portable_setup(self.app)
        ingest_grid(self.app, 0, 9, lambda base, i: base + Decimal("0.4"))
        dispatch(
            self.app,
            {"action": "agreement_run", "specification": base_spec("exp")},
        )

    def test_export_kinds_and_round_trip(self):
        exports = AgreementExports(self.app.agreement)
        pairs_csv = exports.pair_rows("exp")
        self.assertIn("specimen,replicate,reference,candidate", pairs_csv)
        self.assertEqual(pairs_csv.count("\n"), 10)  # header + 9 specimens
        summary_csv = exports.summary_rows("exp")
        self.assertIn("decision", summary_csv)
        exclusions_csv = exports.exclusion_rows("exp")
        self.assertTrue(exclusions_csv.startswith("code,side,observation"))
        payload = json.loads(exports.json_bytes("exp").decode("utf-8"))
        self.assertEqual(payload["specification"]["id"], "exp")

    def test_console_export_writes_files(self):
        with tempfile.TemporaryDirectory() as directory:
            for kind in ("json", "pairs", "exclusions", "summary"):
                destination = str(Path(directory) / f"{kind}.dat")
                result = dispatch(
                    self.app,
                    {
                        "action": "agreement_export",
                        "id": "exp",
                        "kind": kind,
                        "destination": destination,
                    },
                )
                self.assertTrue(Path(destination).is_file())
                self.assertGreater(result["bytes"], 0)
            listing = dispatch(self.app, {"action": "agreement_list"})
            self.assertEqual(
                [item["id"] for item in listing["studies"]], ["exp"]
            )
            fetched = dispatch(
                self.app, {"action": "agreement_get", "id": "exp"}
            )
            self.assertEqual(fetched["id"], "exp")

    def test_unknown_export_kind_rejected(self):
        with self.assertRaises(InputError):
            dispatch(
                self.app,
                {
                    "action": "agreement_export",
                    "id": "exp",
                    "kind": "pdf",
                    "destination": "/tmp/x",
                },
            )


if __name__ == "__main__":
    unittest.main()
