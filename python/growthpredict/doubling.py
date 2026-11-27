"""Successive doubling times of cumulative incidence.

Python port of ``getDoublingTimeCurve.m``.
"""

from __future__ import annotations

import math

import numpy as np

from .errors import GrowthPredictError
from .rt import _colon, _is_real_scalar


def getDoublingTimeCurve(curve, DT, print1=0, priorCumulative=None):
    """Successive doubling times of cumulative incidence.

    Parameters
    ----------
    curve : sequence
        Incidence at successive time points.
    DT : float
        Time between points.
    print1 : 0 or 1
        1 to plot the cumulative curve and the doubling times.
    priorCumulative : float
        Cases counted before the first point of the curve (default 0). For a
        window that starts inside a series, pass the sum of the observations
        before the window.

    Returns
    -------
    tds : ndarray, k-by-2
        Rows ``[time, cumulative count]`` at the start and at each doubling.
    C0data : float or None
        Cumulative count at the first point (None for a curve of one point).
    curve : ndarray
        Cumulative curve at the time points, ``[time, count]``.
    doublingtimes : ndarray
        Time taken by each successive doubling; empty when the count does not
        double within the curve.

    The count that doubles is ``priorCumulative + cumsum(curve)``. Without
    ``priorCumulative`` the count restarts at the first point of the curve, so
    a window taken from the middle of an epidemic shows short 'doublings' of
    its own running total, even when incidence is falling.
    """
    dt1 = 0.01

    if priorCumulative is None or np.size(priorCumulative) == 0:
        priorCumulative = 0

    if not _is_real_scalar(priorCumulative) or not math.isfinite(priorCumulative) or priorCumulative < 0:
        raise GrowthPredictError(
            "GrowthPredict:getDoublingTimeCurve:priorCumulative",
            "priorCumulative must be a number that is zero or positive.",
        )

    tds = np.zeros((0, 2))
    C0data = None
    doublingtimes = np.zeros(0)

    curve = np.asarray(curve, dtype=float).ravel()

    if curve.size <= 1:
        return tds, C0data, curve, doublingtimes

    tb = np.arange(curve.size, dtype=float) * DT

    cumulative1 = priorCumulative + np.cumsum(curve)

    cum1 = cumulative1[-1]

    C0data = float(cumulative1[0])

    XX = _colon(tb[0], dt1, tb[-1])

    cumulative1 = np.interp(XX, tb, cumulative1)

    index1 = cumulative1 <= cum1

    XX = XX[index1]

    cumulative1 = cumulative1[index1]

    step = int(round(1.0 / dt1))
    curveOut = np.column_stack([XX, cumulative1])[::step, :]

    if print1:
        import matplotlib.pyplot as plt

        plt.subplot(1, 2, 1)
        plt.plot(XX, cumulative1, "-")
        plt.xlabel("Time", fontsize=16)
        plt.ylabel("Cumulative cases", fontsize=16)

    x = cumulative1
    timevect = XX

    C0 = x[0]

    rows = [(0.0, float(C0))]

    for _ in range(25):

        C0i1 = 2 * C0

        hits = np.flatnonzero(x <= C0i1)

        if hits.size:
            rows.append((float(timevect[hits[-1]]), float(x[hits[-1]])))

        C0 = C0i1

    tds = np.unique(np.array(rows, dtype=float), axis=0)  # unique(tds,'rows'): sorted rows

    tds = tds[:-1, :]

    if print1:
        import matplotlib.pyplot as plt

        plt.subplot(1, 2, 2)
        if tds.shape[0] > 0:
            plt.plot(tds[:, 1], np.concatenate(([tds[0, 0]], np.diff(tds[:, 0]))), "ro--")
        plt.ylabel("Doubling time", fontsize=16)
        plt.xlabel("Doubling cumulative incidence", fontsize=16)

    if tds.shape[0] <= 1:
        return tds, C0data, curveOut, np.zeros(0)

    maxdoublings = 30

    doublingtimes = np.diff(tds[:, 0])

    doublingtimes = doublingtimes[: min(doublingtimes.size, maxdoublings)]

    return tds, C0data, curveOut, doublingtimes
