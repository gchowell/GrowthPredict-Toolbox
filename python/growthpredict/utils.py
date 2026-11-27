"""Small helpers shared by the toolbox: MATLAB-compatible statistics, the
random-number stream, the data loader and miscellaneous utilities.
"""

from __future__ import annotations

import os

import numpy as np

from .errors import GrowthPredictError

# <============================================================================>
# <============================ Struct container ==============================>
# <============================================================================>


class Struct(dict):
    """A dict whose keys can also be read and written as attributes.

    Stands in for a MATLAB struct: ``s.LB`` and ``s['LB']`` are the same thing.
    """

    __slots__ = ()

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name) from None

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        try:
            del self[name]
        except KeyError:
            raise AttributeError(name) from None

    def copy(self):
        return Struct(self)


# <============================================================================>
# <========================= Random-number stream =============================>
# <============================================================================>

_GLOBAL_RNG = np.random.Generator(np.random.MT19937())


def rng(seed=None, generator="twister"):
    """Seed the toolbox's random-number stream, like MATLAB's ``rng(seed,'twister')``.

    Every function that draws random numbers (start points in ``fit_model``,
    bootstrap noise in ``AddErrorStructure``) uses this stream unless it is
    given its own ``numpy.random.Generator`` through its ``rng`` argument.

    ``rng(1)`` or ``rng(1, 'twister')`` makes a run repeatable in Python. The
    stream is a Mersenne Twister, but NumPy and MATLAB turn its output into
    normal, Poisson and negative binomial variates with different algorithms,
    so a seeded Python run reproduces itself, not the MATLAB run with the same
    seed. Returns the generator.
    """
    global _GLOBAL_RNG
    name = str(generator).lower()
    if name in ("twister", "mt19937", "mt19937ar"):
        bitgen = np.random.MT19937(seed)
    elif name in ("pcg64", "default"):
        bitgen = np.random.PCG64(seed)
    elif name == "philox":
        bitgen = np.random.Philox(seed)
    else:
        raise ValueError("Unknown generator %r; use 'twister', 'pcg64' or 'philox'." % generator)
    _GLOBAL_RNG = np.random.Generator(bitgen)
    return _GLOBAL_RNG


def get_rng(rng_=None):
    """The generator to draw from: the one given, or the toolbox stream."""
    if rng_ is None:
        return _GLOBAL_RNG
    if isinstance(rng_, np.random.Generator):
        return rng_
    return np.random.default_rng(rng_)


# <============================================================================>
# <====================== MATLAB-compatible statistics ========================>
# <============================================================================>


def quantile(x, p, axis=None):
    """Quantiles with MATLAB's definition (``quantile`` / ``prctile``).

    MATLAB treats the sorted sample ``x(1..n)`` as the ``(k-0.5)/n`` quantiles
    and interpolates linearly between them, using the minimum and maximum
    outside that range. This is NumPy's ``method='hazen'``; NumPy's default
    (``'linear'``) gives slightly different interval limits. NaNs are ignored,
    as in MATLAB.

    ``axis=None`` follows MATLAB for a vector (one value per probability). For
    a matrix pass the axis explicitly: MATLAB's ``quantile(X, p)`` works down
    the columns (``axis=0``), ``quantile(X, p, 2)`` along the rows (``axis=1``).
    """
    x = np.asarray(x, dtype=float)
    scalar_p = np.ndim(p) == 0
    p = np.atleast_1d(np.asarray(p, dtype=float))
    if x.size == 0:
        out = np.full(p.shape, np.nan)
        return float(out[0]) if scalar_p else out
    with np.errstate(all="ignore"):
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out = np.nanquantile(x, p, axis=axis, method="hazen")
    if scalar_p:
        out = out[0]
        if np.ndim(out) == 0:
            return float(out)
    return out


def smooth(y, span=5):
    """Moving average like the Curve Fitting Toolbox ``smooth(y, span)``.

    Centred moving average of odd width ``span``; near the ends the window
    shrinks symmetrically so that it stays centred (the first and last values
    are returned unchanged).
    """
    y = np.asarray(y, dtype=float).ravel()
    n = y.size
    if n == 0:
        return y.copy()
    span = int(min(span, n))
    span = span - 1 + span % 2  # force it to be odd
    if span <= 1:
        return y.copy()
    half = (span - 1) // 2
    out = np.empty(n)
    for i in range(n):
        h = min(half, i, n - 1 - i)
        out[i] = y[i - h : i + h + 1].mean()
    return out


def days2weeks(curve):
    """Sum a daily curve into weeks (``days2weeks.m``); the last week may be short."""
    curve = np.asarray(curve, dtype=float).ravel()
    return np.array([curve[i : i + 7].sum() for i in range(0, curve.size, 7)])


def num2str(x):
    """Text of a number as MATLAB's ``num2str`` writes it in the file names."""
    if isinstance(x, (bool, np.bool_)):
        return "1" if x else "0"
    xf = float(x)
    if xf == int(xf):
        return str(int(xf))
    return "%.5g" % xf  # num2str uses about four decimals; adequate for labels


# <============================================================================>
# <=============================== Data loader ================================>
# <============================================================================>

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

# Where the example datasets are looked for when a file is not in ./input:
# a folder 'input' inside the package (present only if the data are bundled
# with an installed copy), then the MATLAB code folder of the repository this
# 'python' folder belongs to.
_BUNDLED_INPUT_DIRS = (
    os.path.join(_PACKAGE_DIR, "input"),
    os.path.join(os.path.dirname(os.path.dirname(_PACKAGE_DIR)), "forecasting_growthmodels code", "input"),
)


def example_input_dir():
    """Folder with the example datasets of the toolbox, or None when it is not found.

    In a checkout of the repository this is ``forecasting_growthmodels code/input``.
    """
    for folder in _BUNDLED_INPUT_DIRS:
        if os.path.isdir(folder):
            return folder
    return None


def resolve_input_file(cadfilename1, input_dir="./input"):
    """Find ``<cadfilename1>.txt``: first in ``input_dir``, then among the example datasets.

    Returns ``(file name with .txt, full path)``. MATLAB looks only in
    ``./input``; the fallback to the example datasets of the repository lets
    the bundled examples run from any working folder.
    """
    name = str(cadfilename1)
    if not name.lower().endswith(".txt"):
        name = name + ".txt"
    candidates = [os.path.join(input_dir, name)] + [os.path.join(d, name) for d in _BUNDLED_INPUT_DIRS]
    for path in candidates:
        if os.path.isfile(path):
            return name, path
    raise FileNotFoundError('File "%s" not found in the specified directory.' % candidates[0])


def load_timeseries(cadfilename1, input_dir="./input"):
    """Load a two-column, headerless text file like MATLAB's ``load``.

    Returns ``(data, file name with .txt)``; ``data`` is an N-by-2 float array
    (time index, incidence). A file whose name begins with ``cumulative``
    (case-insensitive) is differenced into incidence, ``[C(1); diff(C)]``.
    """
    name, path = resolve_input_file(cadfilename1, input_dir)
    data = np.loadtxt(path, ndmin=2)
    if data.size == 0:
        raise GrowthPredictError("GrowthPredict:data", "The dataset is empty")
    if len(name) >= 10 and name[:10].upper() == "CUMULATIVE":
        data = data.copy()
        data[:, 1] = np.concatenate(([data[0, 1]], np.diff(data[:, 1])))  # incidence curve
    return data, name
