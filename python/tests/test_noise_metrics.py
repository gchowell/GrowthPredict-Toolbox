"""Tests of the bootstrap noise, the dispersion estimate, the quantiles, the
forecast scores and the bound diagnostics (no MATLAB test file exists for these).
"""

import math

import numpy as np
import pytest

import growthpredict as gp
from growthpredict import (
    AddErrorStructure,
    GrowthPredictError,
    boundHitFractions,
    computeforecastperformance,
    computeQuantiles,
    computeWIS,
    estimatePearsonDispersion,
    getMeanVarianceRatio,
    quantile,
)

# ---------------------------------------------------------------- quantiles


def test_quantile_follows_matlab_definition():
    # MATLAB documentation example: quantile of x = [2 5 6 10 11 13]
    x = [2, 5, 6, 10, 11, 13]
    assert math.isclose(quantile(x, 0.50), 8.0)
    np.testing.assert_allclose(quantile(x, [0.25, 0.5, 0.75]), [5, 8, 11])
    # sorted values are the (k-0.5)/n quantiles; outside them the extremes are returned
    assert quantile([1, 2, 3, 4], 0.125) == 1
    assert quantile([1, 2, 3, 4], 0.0) == 1
    assert quantile([1, 2, 3, 4], 1.0) == 4
    assert math.isclose(quantile([1, 2, 3, 4], 0.5), 2.5)
    assert math.isclose(quantile([1, 2, 3, 4], 0.25), 1.5)
    assert math.isclose(quantile(np.arange(1, 101), 0.025), 3.0)
    assert math.isclose(quantile(np.arange(1, 101), 0.975), 98.0)
    # NaNs are ignored, as in MATLAB
    assert math.isclose(quantile([1, np.nan, 3], 0.5), 2.0)
    # matrices: along the rows or down the columns
    X = np.arange(12.0).reshape(3, 4)
    np.testing.assert_allclose(quantile(X, 0.5, axis=1), [1.5, 5.5, 9.5])
    np.testing.assert_allclose(quantile(X, 0.5, axis=0), [4, 5, 6, 7])


# ---------------------------------------------------------------- noise

CUM = np.cumsum(np.array([5.0, 20, 60, 150, 300, 420, 380, 200, 90, 30]))
MEANS = np.diff(CUM)


def test_first_point_is_kept_and_shape():
    for dist1 in range(6):
        curves = AddErrorStructure(CUM, 7, dist1, 2.0, 1.5)
        assert curves.shape == (CUM.size, 7)
        assert np.all(curves[0, :] == CUM[0])
        if dist1 > 0:
            assert np.all(curves[1:, :] >= 0)
            assert np.all(curves[1:, :] == np.round(curves[1:, :]))


@pytest.mark.parametrize(
    "dist1,factor1,d,variance",
    [
        (0, 12.0, 0, lambda m: np.full_like(m, 144.0)),
        (1, 1.0, 1, lambda m: m),
        (2, 4.0, 1, lambda m: 4.0 * m),
        (3, 2.5, 1, lambda m: m + 2.5 * m),
        (4, 0.05, 1, lambda m: m + 0.05 * m**2),
        (5, 0.3, 1.5, lambda m: m + 0.3 * m**1.5),
    ],
)
def test_moments_of_the_noise(dist1, factor1, d, variance):
    M = 40000
    curves = AddErrorStructure(CUM, M, dist1, factor1, d, rng=np.random.default_rng(3))
    draws = curves[1:, :]
    target = variance(MEANS)
    se_mean = np.sqrt(target / M)
    assert np.all(np.abs(draws.mean(axis=1) - MEANS) < 5 * se_mean)
    np.testing.assert_allclose(draws.var(axis=1, ddof=1), target, rtol=0.08)


def test_negative_binomial_without_overdispersion_is_poisson_and_zero_mean_is_handled():
    flat = np.array([10.0, 10, 10, 60])  # two zero means
    curves = AddErrorStructure(flat, 20000, 2, 1.0, 1, rng=np.random.default_rng(1))
    assert np.all(np.isfinite(curves))
    assert curves[1, :].mean() < 0.01  # mean replaced by 0.001
    np.testing.assert_allclose(curves[3, :].var(ddof=1), 50, rtol=0.05)  # VAR = mean: Poisson


def test_noise_errors():
    with pytest.raises(GrowthPredictError) as e:
        AddErrorStructure(CUM, 2, 6, 1, 1)
    assert e.value.identifier == "GrowthPredict:AddErrorStructure:UnsupportedDist"
    with pytest.raises(GrowthPredictError) as e:
        AddErrorStructure(CUM, 2, 0, float("nan"), 0)
    assert e.value.identifier == "GrowthPredict:AddErrorStructure:InvalidDraws"


def test_seed_makes_draws_repeatable():
    gp.rng(5)
    a = AddErrorStructure(CUM, 3, 1, 1, 1)
    gp.rng(5)
    b = AddErrorStructure(CUM, 3, 1, 1, 1)
    c = AddErrorStructure(CUM, 3, 1, 1, 1)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


def test_pearson_dispersion():
    mu = np.array([10.0, 20, 40, 80, 0.2, 50])
    y = np.array([10.0, 26, 34, 92, 3, 41])
    # fixI0 = 1 drops the first point, the fitted mean 0.2 is below 0.5
    use = [1, 2, 3, 5]
    expected = float(np.sum((y[use] - mu[use]) ** 2 / mu[use]) / (4 - 2))
    phi, info = estimatePearsonDispersion(y, mu, 2, 1, 1)
    assert math.isclose(phi, expected)
    assert (info.nUsed, info.dof, info.excludedSmallMean, info.excludedFixed) == (4, 2, 1, 1)

    # too few usable points
    phi, info = estimatePearsonDispersion(y, mu, 4, 1, 1)
    assert math.isnan(phi) and info.dof == 0

    # Poisson data around the true curve: dispersion close to 1
    rng = np.random.default_rng(0)
    mu = 200 * np.exp(-0.5 * ((np.arange(400) - 200) / 60.0) ** 2) + 5
    phi, _ = estimatePearsonDispersion(rng.poisson(mu), mu, 3)
    assert abs(phi - 1) < 0.2

    with pytest.raises(GrowthPredictError):
        estimatePearsonDispersion([1, 2, 3], [1, 2], 1)
    with pytest.raises(GrowthPredictError):
        estimatePearsonDispersion([1, 2, 3], [1, 2, 3], 1, 2)


def test_mean_variance_ratio():
    x = np.array([1.0, 3, 5, 2, 4, 9, 7, 7])
    ratios, bins = getMeanVarianceRatio(x, 4, 1)
    assert ratios.shape == (2, 3)
    assert math.isclose(ratios[0, 1], 2.75) and math.isclose(ratios[0, 2], np.var(x[:4], ddof=1))
    assert math.isclose(ratios[0, 0], ratios[0, 2] / ratios[0, 1])
    assert bins.tolist() == [2.5, 6.5]
    ratios, bins = getMeanVarianceRatio(x, 4, 2)
    assert ratios.shape == (5, 3) and bins.tolist() == [1, 2, 3, 4, 5]


# ---------------------------------------------------------------- scores


def reference_wis(y, sample):
    """Weighted interval score written out as in Bracher et al. (2021), one observation."""
    alphas = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    total = 0.5 * abs(y - max(quantile(sample, 0.5), -np.inf))
    for alpha in alphas:
        lo = max(quantile(sample, alpha / 2), 0)
        hi = max(quantile(sample, 1 - alpha / 2), 0)
        score = (hi - lo) + (2 / alpha) * max(lo - y, 0) + (2 / alpha) * max(y - hi, 0)
        total += (alpha / 2) * score
    return total / (len(alphas) + 0.5)


def test_wis_coverage_and_point_errors():
    rng = np.random.default_rng(11)
    n, h, ncol = 12, 4, 400
    truth = 100 + 10 * np.arange(n + h)
    noisy = truth[:, None] + rng.normal(0, 15, size=(n + h, ncol))
    smooth = truth[:, None] + rng.normal(0, 2, size=(n + h, 50))
    obs = truth + rng.normal(0, 15, size=n + h)
    obs[n + 1] = truth[n + 1] + 200  # one forecast point far outside the interval
    datalatest = np.column_stack([np.arange(n + h), obs])
    data1 = datalatest[:n, :]

    WISC, WISFS, qc, qf = computeWIS(data1, datalatest, noisy, h)

    expected_c = np.mean([reference_wis(obs[i], noisy[i, :]) for i in range(n)])
    assert math.isclose(WISC, expected_c, rel_tol=1e-12)
    wf = [reference_wis(obs[n + i], noisy[n + i, :]) for i in range(h)]
    assert WISFS.shape == (h, 2)
    np.testing.assert_allclose(WISFS[:, 0], [1, 2, 3, 4])
    np.testing.assert_allclose(WISFS[:, 1], np.cumsum(wf) / np.arange(1, h + 1), rtol=1e-12)
    assert qc.shape == (n, 23) and qf.shape == (h, 23)

    out = computeforecastperformance(data1, datalatest, smooth, noisy, h)
    RMSECS, MSECS, MAECS, PICS, MISCS, RMSEFS, MSEFS, MAEFS, PIFS, MISFS = out
    med = np.median(smooth, axis=1)
    assert math.isclose(MAECS[-1, -1], np.mean(np.abs(obs[:n] - med[:n])))
    assert math.isclose(MSECS[-1, -1], np.mean((obs[:n] - med[:n]) ** 2))
    assert math.isclose(RMSECS[-1, -1], math.sqrt(MSECS[-1, -1]))
    # horizon h summarizes ALL of the first h held-out observations
    for k in range(1, h + 1):
        assert math.isclose(MAEFS[k - 1, 1], np.mean(np.abs(obs[n : n + k] - med[n : n + k])))
        assert math.isclose(MSEFS[k - 1, 1], np.mean((obs[n : n + k] - med[n : n + k]) ** 2))
    lo = np.array([quantile(noisy[i, :], 0.025) for i in range(n + h)])
    hi = np.array([quantile(noisy[i, :], 0.975) for i in range(n + h)])
    inside = (obs >= lo) & (obs <= hi)
    assert math.isclose(PICS[-1, -1], 100 * inside[:n].mean())
    np.testing.assert_allclose(PIFS[:, 1], 100 * np.cumsum(inside[n:]) / np.arange(1, h + 1))
    assert PIFS[1, 1] == 50.0 and PIFS[0, 1] == 100.0  # the outlier is the second forecast point
    mis = (hi - lo) + 40 * (lo - obs) * (obs < lo) + 40 * (obs - hi) * (obs > hi)
    assert math.isclose(MISCS[-1, -1], mis[:n].mean())
    np.testing.assert_allclose(MISFS[:, 1], np.cumsum(mis[n:]) / np.arange(1, h + 1))

    quantilesc, quantilesf, alphaquantiles = computeQuantiles(data1, noisy, h)
    assert alphaquantiles.size == 23 and alphaquantiles[0] == 0.01 and alphaquantiles[-1] == 0.99
    np.testing.assert_allclose(quantilesc, qc)
    np.testing.assert_allclose(quantilesf, qf)
    np.testing.assert_allclose(quantilesc[:, 1], lo[:n])
    np.testing.assert_allclose(quantilesf[:, 21], hi[n:])
    assert np.all(np.diff(quantilesc, axis=1) >= 0)


def test_scores_without_forecast_or_without_data():
    rng = np.random.default_rng(2)
    curves = np.abs(rng.normal(50, 5, size=(8, 100)))
    data = np.column_stack([np.arange(6), np.full(6, 50.0)])

    WISC, WISFS, qc, qf = computeWIS(data, data, curves[:6], 0)
    assert WISFS.shape == (0, 2) and qf.shape == (0, 23)
    out = computeforecastperformance(data, data, curves[:6], curves[:6], 0)
    assert all(a.shape == (0, 2) for a in out[5:])

    # forecast rows exist but the observations do not: calibration metrics only
    with pytest.warns(UserWarning, match="too long"):
        WISC2, WISFS, _, qf = computeWIS(data, data, curves, 2)
    assert math.isclose(WISC2, WISC) and WISFS.shape == (0, 2) and qf.shape == (2, 23)
    with pytest.warns(UserWarning, match="too long"):
        out = computeforecastperformance(data, data, curves, curves, 2)
    assert all(a.shape == (0, 2) for a in out[5:])

    # negative prediction limits are set to zero
    neg = rng.normal(0, 5, size=(6, 200))
    qc, _, _ = computeQuantiles(data, neg, 0)
    assert qc.min() == 0


def test_bound_hit_fractions():
    P = np.array([[1.0, 5.0, 2.0], [1.0005, 9.0, 2.0], [3.0, 9.999, 2.0], [np.nan, 7.0, 2.0]])
    LB = [1.0, 0.0, 2.0]
    UB = [10.0, 10.0, 2.0]
    fracLB, fracUB, hitLB, hitUB = boundHitFractions(P, LB, UB)
    assert math.isclose(fracLB[0], 2 / 3)  # NaN draws are left out; 1.0005 is within 0.1% of LB
    assert fracUB[0] == 0
    assert fracLB[1] == 0 and math.isclose(fracUB[1], 1 / 4)  # 9.999 is within 0.1% of UB
    assert math.isnan(fracLB[2]) and math.isnan(fracUB[2])  # pinned column
    assert hitLB.shape == P.shape and hitUB.shape == P.shape

    # a box longer than the parameters, per-draw boxes and a custom tolerance
    f1, f2, _, _ = boundHitFractions(P[:, :2], LB + [0, 0], UB + [1, 1])
    np.testing.assert_allclose(f1, fracLB[:2])
    f1, f2, _, _ = boundHitFractions(P[:, :2], np.tile(LB[:2], (4, 1)), np.tile(UB[:2], (4, 1)), 0.2)
    assert math.isclose(f2[1], 2 / 4)

    # a bound of exactly zero counts exact hits only
    f1, _, _, _ = boundHitFractions(np.array([[0.0], [1e-9], [0.5]]), [0.0], [1.0])
    assert math.isclose(f1[0], 1 / 3)

    with pytest.raises(GrowthPredictError) as e:
        boundHitFractions(P, [1.0, 0.0], UB)
    assert e.value.identifier == "boundHitFractions:sizeMismatch"
