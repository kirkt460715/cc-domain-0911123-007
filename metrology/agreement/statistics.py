"""Method agreement statistics.

All numeric work uses :class:`decimal.Decimal` so results are reproducible
without third-party packages.  The module is deliberately free of database or
ledger concerns; :mod:`metrology.agreement.studies` supplies the paired
observations and interprets these results against the quality requirements.
"""

import math
import random
from dataclasses import dataclass, field
from decimal import Decimal, localcontext

from ..errors import InputError
from ..values import decimal_text, number

D = Decimal
ZERO = D(0)
ONE = D(1)
TWO = D(2)
Z95 = D("1.95996398454005423552")
T95 = {  # two-sided 95% Student-t critical values, degrees of freedom -> t
    1: D("12.706204736432095"),
    2: D("4.302652729911275"),
    3: D("3.182446305281705"),
    4: D("2.7764451051977987"),
    5: D("2.5705818366147395"),
    6: D("2.4469118511449717"),
    7: D("2.3646242515927846"),
    8: D("2.3060041350333754"),
    9: D("2.262157162740992"),
    10: D("2.2281388519649385"),
    15: D("2.131449545559323"),
    20: D("2.085963447265837"),
    25: D("2.059538552753296"),
    30: D("2.042272455498902"),
    40: D("2.021075382995467"),
    60: D("2.0002978210582616"),
    120: D("1.979930400814168"),
}
SCALES = ("absolute", "percent", "log")
MAX_BOOTSTRAP = 20000
MAX_POINTS = 10000


def sqrt(value):
    if value == 0:
        return ZERO
    with localcontext() as ctx:
        ctx.prec = 40
        return value.sqrt()


def mean(values):
    values = list(values)
    if not values:
        return None
    return sum(values, ZERO) / D(len(values))


def variance(values, mu=None):
    """Population variance (divide by n); empty or singleton gives zero."""
    values = list(values)
    if len(values) < 2:
        return ZERO
    mu = mean(values) if mu is None else mu
    return sum(((v - mu) ** 2 for v in values), ZERO) / D(len(values))


def t_critical(degrees):
    if degrees in T95:
        return T95[degrees]
    for key in sorted(T95):
        if key >= degrees:
            return T95[key]
    return Z95


@dataclass(frozen=True)
class Point:
    """One paired reading in the analysis unit of the study."""

    specimen: str
    replicate: str
    reference: Decimal
    candidate: Decimal
    reference_id: str
    candidate_id: str
    reference_revision: int
    candidate_revision: int
    certificate_reference: str
    certificate_candidate: str
    reference_at: str
    candidate_at: str

    @property
    def average(self):
        return (self.reference + self.candidate) / TWO

    @property
    def difference(self):
        return self.candidate - self.reference


@dataclass(frozen=True)
class Transformed:
    scale: str
    points: tuple
    # transformed reference/candidate per point, and transformed difference
    references: tuple
    candidates: tuple
    differences: tuple
    averages: tuple
    diagnostics: tuple = field(default_factory=tuple)


def _ln(value):
    with localcontext() as ctx:
        ctx.prec = 40
        return D(repr(math.log(float(value))))


def transform(points, scale, percent_fraction=None):
    scale = scale.lower()
    if scale not in SCALES:
        raise InputError("unsupported difference scale", field="scale")
    diagnostics = []
    references, candidates, differences, averages = [], [], [], []
    for point in points:
        ref, cand = point.reference, point.candidate
        if scale == "absolute":
            diff = cand - ref
            avg = (cand + ref) / TWO
        elif scale == "percent":
            fraction = number(percent_fraction, "percent_fraction", positive=True)
            if ref <= 0:
                raise InputError(
                    "percent scale requires strictly positive reference values",
                    field="scale",
                )
            diff = (cand - ref) / ref * fraction
            avg = (cand + ref) / TWO
        else:  # log
            if ref <= 0 or cand <= 0:
                raise InputError(
                    "log scale requires strictly positive values on both methods",
                    field="scale",
                )
            ref, cand = _ln(ref), _ln(cand)
            diff = cand - ref
            avg = (cand + ref) / TWO
        references.append(ref)
        candidates.append(cand)
        differences.append(diff)
        averages.append(avg)
    if scale == "log":
        diagnostics.append(
            {
                "code": "log_scale",
                "message": "differences are natural-log ratios; exp(limits) gives "
                "multiplicative ratio bounds",
            }
        )
    return Transformed(
        scale,
        tuple(points),
        tuple(references),
        tuple(candidates),
        tuple(differences),
        tuple(averages),
        tuple(diagnostics),
    )


def pearson(xs, ys):
    if len(xs) < 2:
        return None
    mx, my = mean(xs), mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / sqrt(sxx * syy)


def _specimen_groups(points):
    groups = {}
    for point in points:
        groups.setdefault(point.specimen, []).append(point)
    return {key: tuple(value) for key, value in sorted(groups.items())}


def bland_altman(transformed):
    """Bland-Altman summary.

    With repeated measurements the standard deviation of the differences is
    split into between-specimen and within-specimen components using a
    one-way random-effects identity on the differences; the repeatability
    component is the pooled within-specimen variance.  Limits use the total
    variance so they describe agreement of individual readings, while the
    repeatability limit uses only the within-specimen component.
    """
    points = transformed.points
    differences = list(transformed.differences)
    if not points:
        return None
    groups = _specimen_groups(points)
    n = len(points)
    bias = mean(differences)
    total_variance = variance(differences, bias)
    within_terms = []
    within_df = 0
    specimen_means = []
    transformed_by_id = {
        id(point): diff for point, diff in zip(points, differences)
    }
    for group in groups.values():
        group_differences = [transformed_by_id[id(p)] for p in group]
        group_mu = mean(group_differences)
        specimen_means.append(group_mu)
        if len(group) > 1:
            within_terms.append(
                sum((d - group_mu) ** 2 for d in group_differences)
            )
            within_df += len(group) - 1
    within_variance = (
        sum(within_terms, ZERO) / D(within_df) if within_df else ZERO
    )
    between_variance_raw = total_variance - within_variance
    between_variance = max(between_variance_raw, ZERO)
    repeats = n > len(groups)
    if repeats and within_df > 0:
        critical = t_critical(within_df)
        repeatability_limit = critical * sqrt(TWO * within_variance)
    else:
        repeatability_limit = None
    sd_total = sqrt(total_variance)
    lower = bias - Z95 * sd_total
    upper = bias + Z95 * sd_total
    result = {
        "pairs": n,
        "specimens": len(groups),
        "repeated_specimens": sum(1 for group in groups.values() if len(group) > 1),
        "scale": transformed.scale,
        "bias": decimal_text(bias),
        "sd_difference": decimal_text(sd_total),
        "lower_limit": decimal_text(lower),
        "upper_limit": decimal_text(upper),
        "mean_of_averages": decimal_text(mean(transformed.averages)),
        "correlation": None,
        "components": {
            "between_specimen_variance": decimal_text(between_variance),
            "within_specimen_variance": decimal_text(within_variance),
            "negative_variance_clipped": between_variance_raw < 0,
        },
        "repeatability_limit": (
            decimal_text(repeatability_limit)
            if repeatability_limit is not None
            else None
        ),
        "repeatability_df": within_df,
    }
    correlation = pearson(transformed.references, transformed.candidates)
    result["correlation"] = (
        decimal_text(correlation) if correlation is not None else None
    )
    # proportional-bias regression of difference on average (descriptive only).
    avgs = list(transformed.averages)
    ma = mean(avgs)
    sxx = sum((a - ma) ** 2 for a in avgs)
    if sxx > 0:
        slope = sum(
            (a - ma) * (d - bias) for a, d in zip(avgs, differences)
        ) / sxx
        intercept = bias - slope * ma
        result["proportional_bias"] = {
            "slope": decimal_text(slope),
            "intercept": decimal_text(intercept),
        }
    else:
        result["proportional_bias"] = None
    return result


def deming(references, candidates, variance_ratio):
    """Deming regression candidate = slope * reference + intercept.

    ``variariance_ratio`` is the candidate-method error variance divided by
    the reference-method error variance (lambda).  Ordinary least squares is
    never used: lambda must be a finite positive number.
    """
    ratio = number(variance_ratio, "variance_ratio", positive=True)
    refs = [number(v) for v in references]
    cands = [number(v) for v in candidates]
    if len(refs) < 3:
        raise InputError("Deming regression requires at least three specimens")
    if len(set(refs)) < 2 or len(set(cands)) < 2:
        raise InputError("degenerate data cannot be regressed")
    mx, my = mean(refs), mean(cands)
    sxx = sum((x - mx) ** 2 for x in refs) / D(len(refs))
    syy = sum((y - my) ** 2 for y in cands) / D(len(refs))
    sxy = sum((x - mx) * (y - my) for x, y in zip(refs, cands)) / D(len(refs))
    if sxx <= 0 or sxy == 0:
        raise InputError("singular Deming regression")
    discriminant = (syy - ratio * sxx) ** 2 + FOUR * sxy * sxy
    slope = (syy - ratio * sxx + sqrt(discriminant)) / (TWO * sxy)
    intercept = my - slope * mx
    # Jackknife intervals over the supplied points.
    slopes, intercepts = [], []
    for index in range(len(refs)):
        sub_x = refs[:index] + refs[index + 1 :]
        sub_y = cands[:index] + cands[index + 1 :]
        try:
            s, i = _deming_coefficients(sub_x, sub_y, ratio)
        except InputError:
            continue
        slopes.append(s)
        intercepts.append(i)
    intervals = {}
    if len(slopes) >= 3:
        intervals = {
            "slope_interval": _jackknife_interval(slopes),
            "intercept_interval": _jackknife_interval(intercepts),
        }
    return {
        "variance_ratio": decimal_text(ratio),
        "slope": decimal_text(slope),
        "intercept": decimal_text(intercept),
        "n": len(refs),
        **intervals,
    }


FOUR = D(4)


def _deming_coefficients(refs, cands, ratio):
    mx, my = mean(refs), mean(cands)
    sxx = sum((x - mx) ** 2 for x in refs) / D(len(refs))
    syy = sum((y - my) ** 2 for y in cands) / D(len(refs))
    sxy = sum((x - mx) * (y - my) for x, y in zip(refs, cands)) / D(len(refs))
    if sxx <= 0 or sxy == 0:
        raise InputError("singular Deming regression")
    discriminant = (syy - ratio * sxx) ** 2 + FOUR * sxy * sxy
    slope = (syy - ratio * sxx + sqrt(discriminant)) / (TWO * sxy)
    return slope, my - slope * mx


def _jackknife_interval(values):
    n = D(len(values))
    mu = mean(values)
    if n < 2:
        return None
    se = sqrt((n - ONE) / n * sum((v - mu) ** 2 for v in values))
    critical = t_critical(len(values) - 1)
    return {
        "lower": decimal_text(mu - critical * se),
        "upper": decimal_text(mu + critical * se),
        "standard_error": decimal_text(se),
    }


def cluster_bootstrap(transformed, statistic, *, replicates, seed):
    """Resample whole specimens (all replicates travel together).

    ``statistic`` accepts a list of :class:`Point` and returns a dict of
    named Decimal scalars.  Returns percentile intervals per statistic.
    """
    replicates = int(replicates)
    if not 10 <= replicates <= MAX_BOOTSTRAP:
        raise InputError("bootstrap replicates must be between ten and twenty thousand")
    groups = list(_specimen_groups(transformed.points).values())
    if len(groups) < 2:
        raise InputError("cluster bootstrap needs at least two specimens")
    rng = random.Random(seed)
    collected = {}
    for _ in range(replicates):
        sampled = []
        for index in (rng.randrange(len(groups)) for _ in groups):
            sampled.extend(groups[index])
        try:
            values = statistic(sampled)
        except InputError:
            continue
        for key, value in values.items():
            if value is None:
                continue
            collected.setdefault(key, []).append(number(value))
    return {
        key: {
            "lower": decimal_text(_percentile(entries, D("0.025"))),
            "upper": decimal_text(_percentile(entries, D("0.975"))),
            "replicates": len(entries),
        }
        for key, entries in sorted(collected.items())
    }


def _percentile(values, probability):
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    low = int(position)
    high = min(low + 1, len(ordered) - 1)
    fraction = position - low
    return ordered[low] + fraction * (ordered[high] - ordered[low])


def ba_bootstrap_statistic(scale, percent_fraction):
    """Statistic callback producing BA scalars on a resampled point list."""

    def statistic(points):
        transformed = transform(points, scale, percent_fraction)
        summary = bland_altman(transformed)
        return {
            "bias": number(summary["bias"]),
            "lower_limit": number(summary["lower_limit"]),
            "upper_limit": number(summary["upper_limit"]),
            "sd_difference": number(summary["sd_difference"]),
        }

    return statistic
