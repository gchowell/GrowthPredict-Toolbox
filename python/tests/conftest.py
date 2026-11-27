"""Shared fixtures of the GrowthPredict test suite.

Run from the ``python`` folder::

    python -m pytest            # everything (a few minutes)
    python -m pytest -m "not slow"   # without the end-to-end runs
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import growthpredict as gp  # noqa: E402


@pytest.fixture(scope="session")
def mpox():
    """The bundled U.S. mpox series: N-by-2 array [time index, weekly cases]."""
    data, _ = gp.load_timeseries("Most_Recent_Timeseries_US-CDC")
    return data


@pytest.fixture(autouse=True)
def _seed_and_close():
    """Every test starts from rng(1) and leaves no figure open."""
    import matplotlib.pyplot as plt

    gp.rng(1)
    yield
    plt.close("all")


def error_id(excinfo):
    """Identifier of a GrowthPredictError caught by pytest.raises."""
    return excinfo.value.identifier


def assert_rel(actual, expected, rtol):
    np.testing.assert_allclose(actual, expected, rtol=rtol, atol=0)
