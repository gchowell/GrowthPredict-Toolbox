"""Regression tests for the time origin of the Gompertz model (port of tests/test_gompertz.m).

The Gompertz model dC/dt = r*C*exp(-a*t) is the only model whose growth rate
depends on time. Time is measured from the first point of the curve, so the
result does not depend on the time label of the first observation.

Reference values are those of the MATLAB test file (least squares).
"""

import numpy as np

import growthpredict as gp
from growthpredict import fit_model, growthIncidenceCurve, initialParams


def local_fit(data1, DT):
    params0 = initialParams(data1[:, 1], 5)
    res = fit_model(data1, params0, 1, 10, DT, 5, 0, None, method1=0)
    return res.P, res.fval, res.fitDiagnostics


def test_curve_matches_closed_form():
    # C(tau) = C0*exp((r/a)*(1 - exp(-a*tau))), tau = time since the first point.
    r, a, C0 = 0.8, 0.15, 20

    tau = np.arange(10.0)
    C = C0 * np.exp((r / a) * (1 - np.exp(-a * tau)))
    expected = np.concatenate(([C[0]], np.diff(C)))

    incidence = growthIncidenceCurve([r, 1, a, 1, C0], tau, 5)

    np.testing.assert_allclose(incidence, expected, atol=0.05)
    assert abs(incidence[1] - 22.0401) < 0.01  # reference

    # the port solves much more tightly than that
    np.testing.assert_allclose(incidence, expected, rtol=1e-6)


def test_curve_does_not_depend_on_time_labels():
    # The same parameters give the same curve for labels 0..19, 12..31 and 100..119.
    z = [0.8, 1, 0.15, 1, 20]

    base = growthIncidenceCurve(z, np.arange(20.0), 5)

    for shift1 in (1, 12, 100):
        shifted = growthIncidenceCurve(z, np.arange(20.0) + shift1, 5)
        np.testing.assert_allclose(shifted, base, rtol=1e-9)

    # The other models never depended on the labels.
    zGLM = [0.9, 0.8, 1, 5000, 20]

    base = growthIncidenceCurve(zGLM, np.arange(20.0), 1)
    shifted = growthIncidenceCurve(zGLM, np.arange(12.0, 32.0), 1)

    np.testing.assert_allclose(shifted, base, atol=0.05)


def test_fit_does_not_depend_on_time_labels(mpox):
    # U.S. mpox, rows 1-20, labelled 0..19, 12..31 and 100..119.
    y = mpox[:20, 1]

    for shift1 in (0, 12, 100):

        gp.rng(1)

        P, fval, diagnostics = local_fit(np.column_stack([np.arange(20.0) + shift1, y]), 1)

        np.testing.assert_allclose(P[0], 1.4738, rtol=0.02)  # reference
        np.testing.assert_allclose(P[2], 0.18763, rtol=0.02)  # reference
        np.testing.assert_allclose(fval, 2.71708e6, rtol=0.01)  # reference
        assert not (diagnostics.atLB[0] or diagnostics.atUB[0])


def test_windows_inside_the_series(mpox):
    # Rows 8-27 (labels 7..26) and rows 13-32 (labels 12..31).
    P, fval, diagnostics = local_fit(mpox[7:27, :], 1)

    np.testing.assert_allclose(P[0], 0.99825, rtol=0.02)  # reference
    np.testing.assert_allclose(P[2], 0.30114, rtol=0.02)  # reference
    np.testing.assert_allclose(fval, 185716, rtol=0.01)  # reference
    assert not (diagnostics.atLB[0] or diagnostics.atUB[0])

    P, fval, diagnostics = local_fit(mpox[12:32, :], 1)

    np.testing.assert_allclose(P[0], 0.6732, rtol=0.02)  # reference
    np.testing.assert_allclose(P[2], 0.3939, rtol=0.02)  # reference
    np.testing.assert_allclose(fval, 452472, rtol=0.01)  # reference
    assert not (diagnostics.atLB[0] or diagnostics.atUB[0])


def test_time_step_rescales_the_rates(mpox):
    # With DT = 7 the same observations give the same curve and rates 7 times
    # smaller. The bounds of r and a are per unit of time.
    data1 = np.column_stack([np.arange(20.0), mpox[:20, 1]])

    P, fval, diagnostics = local_fit(data1, 7)

    np.testing.assert_allclose(7 * P[0], 1.4738, rtol=0.02)
    np.testing.assert_allclose(7 * P[2], 0.18763, rtol=0.02)
    np.testing.assert_allclose(fval, 2.71708e6, rtol=0.01)
    np.testing.assert_allclose(7 * diagnostics.UB[0], 10.005, rtol=1e-3)  # reference


def test_forecast_continues_the_fit(mpox):
    # The forecast curve is solved from the same origin as the fitted curve.
    data1 = mpox[12:32, :]

    params0 = initialParams(data1[:, 1], 5)

    res = fit_model(data1, params0, 1, 10, 1, 5, 4, None, method1=0)

    assert res.forecastcurve.size == 24
    assert res.timevect2.size == 24
    np.testing.assert_allclose(res.forecastcurve[:20], res.fitcurve, atol=0.5)
    np.testing.assert_allclose(res.forecastcurve[20:24], [5.22, 3.52, 2.38, 1.60], atol=0.3)  # reference
