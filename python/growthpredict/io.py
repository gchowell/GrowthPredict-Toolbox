"""Output files: names, MAT files and CSV tables.

The file names are those of the MATLAB toolbox, so results written by the
Python port sit next to (and are named like) results written by MATLAB.

MAT files are written with ``scipy.io.savemat`` (MATLAB v5 format, readable
with MATLAB's ``load``) and hold the variables listed in the README of the
toolbox under the same names (``P_model1d``, ``Phatss_model1``,
``fitDiagnostics_model1d``, ``forecast_model12``, ``timevect2`` ...). CSV
tables have the column names of MATLAB's ``writetable`` output and numbers
with 15 significant digits.
"""

from __future__ import annotations

import os

import numpy as np
from scipy.io import loadmat, savemat

from .metrics import ALPHAQUANTILES
from .utils import Struct, num2str

QUANT_NAMES_RANKED = ["Q_%.3f" % q for q in ALPHAQUANTILES]


class OutputNames:
    """Builds the output file names of one analysis (same pattern as MATLAB)."""

    def __init__(
        self, output_dir, cadfilename1, flag1, fixI0, method1, dist1, windowsize1, forecastingperiod, caddisease, datatype
    ):
        self.output_dir = output_dir
        self.cadfilename1 = cadfilename1  # with '.txt', as in MATLAB
        self.flag1 = flag1
        self.fixI0 = fixI0
        self.method1 = method1
        self.dist1 = dist1
        self.windowsize1 = windowsize1
        self.forecastingperiod = forecastingperiod
        self.caddisease = caddisease
        self.datatype = datatype

    def _path(self, name):
        return os.path.join(self.output_dir, name)

    def _model(self):
        return "-flag1-%s-fixI0-%s-method-%s-dist-%s" % (
            num2str(self.flag1), num2str(self.fixI0), num2str(self.method1), num2str(self.dist1),
        )  # fmt: skip

    def _tail(self):
        return "-calibrationperiod-%s-horizon-%s-%s-%s.csv" % (
            num2str(self.windowsize1), num2str(self.forecastingperiod), self.caddisease, self.datatype,
        )  # fmt: skip

    # --- MAT files -----------------------------------------------------------

    def _mat(self, prefix, tstart, tend, forecastingperiod=None):
        fp = self.forecastingperiod if forecastingperiod is None else forecastingperiod
        return self._path(
            "%s-growthModel-%s%s-tstart-%s-tend-%s-calibrationperiod-%s-forecastingperiod-%s.mat"
            % (prefix, self.cadfilename1, self._model(), num2str(tstart), num2str(tend), num2str(self.windowsize1), num2str(fp))
        )

    def window_mat(self, i, tend1, forecastingperiod=None):
        """Per-window results: ``Forecast-growthModel-<file>-...-tstart-<i>-tend-<tend1>-....mat``."""
        return self._mat("Forecast", i, tend1, forecastingperiod)

    def summary_mat(self, prefix, tstart1, tend1):
        """Rolling summaries: ``parameters-``, ``performanceCalibration-``, ``QuantilesCalibration-`` ..."""
        return self._mat(prefix, tstart1, tend1)

    # --- CSV files -----------------------------------------------------------

    def rolling_csv(self, prefix, tstart1, tend1):
        """``<prefix>-flag1-..-tstart-..-tend-..-calibrationperiod-..-horizon-..-<disease>-<datatype>.csv``."""
        return self._path(
            "%s%s-tstart-%s-tend-%s%s" % (prefix, self._model(), num2str(tstart1), num2str(tend1), self._tail())
        )

    def window_csv_tstart_first(self, prefix, i):
        """``<prefix>-flag1-<f>-tstart-<i>-fixI0-..-method-..-dist-..-calibrationperiod-...csv``
        (Fit-, quantile-fit-, quantile-, doublingtimes-, Rt- files)."""
        return self._path(
            "%s-flag1-%s-tstart-%s-fixI0-%s-method-%s-dist-%s%s"
            % (prefix, num2str(self.flag1), num2str(i), num2str(self.fixI0), num2str(self.method1),
               num2str(self.dist1), self._tail())
        )  # fmt: skip

    def window_csv_tstart_last(self, prefix, i):
        """``<prefix>-flag1-..-fixI0-..-method-..-dist-..-tstart-<i>-calibrationperiod-...csv``
        (Forecast-, quantile-forecast- files)."""
        return self._path("%s%s-tstart-%s%s" % (prefix, self._model(), num2str(i), self._tail()))

    def rts_csv(self, tstart1, tend1):
        """``Rts-flag1-<f>-tstart-..-tend-..-fixI0-..-method-..-dist-..-calibrationperiod-...csv`` (plotRt)."""
        return self._path(
            "Rts-flag1-%s-tstart-%s-tend-%s-fixI0-%s-method-%s-dist-%s%s"
            % (num2str(self.flag1), num2str(tstart1), num2str(tend1), num2str(self.fixI0),
               num2str(self.method1), num2str(self.dist1), self._tail())
        )  # fmt: skip


# <============================================================================>
# <================================ MAT files =================================>
# <============================================================================>


def _matlab_value(v):
    """Make a Python value storable by savemat (None -> [], nested dicts kept)."""
    if v is None:
        return np.zeros((0, 0))
    if isinstance(v, dict):
        return {k: _matlab_value(x) for k, x in v.items()}
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (list, tuple)):
        try:
            return np.asarray(v, dtype=float)
        except (TypeError, ValueError):
            return np.asarray(v, dtype=object)
    return v


def save_mat(path, variables):
    """Write a dict of variables to a MAT file (parameter vectors as rows, like MATLAB)."""
    out = {}
    for key, value in variables.items():
        out[key] = _matlab_value(value)
    savemat(path, out, do_compression=True, oned_as="row", long_field_names=True)


def load_mat(path):
    """Read a MAT file written by :func:`save_mat` into a :class:`Struct`.

    Numeric arrays come back two-dimensional, as MATLAB stores them; use
    :func:`scalar` and :func:`vector` to unwrap them.
    """
    raw = loadmat(path, squeeze_me=False, struct_as_record=False, mat_dtype=True)
    return Struct({k: v for k, v in raw.items() if not k.startswith("__")})


def scalar(x):
    """The single number held by a 1-by-1 MAT array."""
    return float(np.asarray(x, dtype=float).ravel()[0])


def vector(x):
    """A MAT row or column as a one-dimensional float array."""
    return np.asarray(x, dtype=float).ravel()


# <============================================================================>
# <================================ CSV tables ================================>
# <============================================================================>


def _fmt(x):
    x = float(x)
    if np.isnan(x):
        return "NaN"
    if np.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return "%.15g" % x


def writetable(path, array, names):
    """Write a numeric table with a header row, like ``writetable(array2table(A), file)``."""
    array = np.asarray(array, dtype=float)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if array.size and array.shape[1] != len(names):
        raise ValueError("The table has %d columns but %d names were given." % (array.shape[1], len(names)))
    with open(path, "w", newline="") as fh:
        fh.write(",".join(names) + "\n")
        for row in array:
            fh.write(",".join(_fmt(v) for v in row) + "\n")


def readtable(path):
    """Read a CSV written by :func:`writetable`. Returns ``(names, array)``."""
    with open(path, "r") as fh:
        names = fh.readline().rstrip("\r\n").split(",")
    array = np.genfromtxt(path, delimiter=",", skip_header=1, ndmin=2)
    return names, array


def _triples(label):
    return ["%s estimate" % label, "%s LB" % label, "%s UB" % label]


def write_rolling_tables(names, tstart1, tend1, method1, s):
    """``parameters-``, ``MCSES-`` and ``SCIS-rollingwindow`` CSV files."""
    time = np.arange(tstart1, tend1 + 1, dtype=float).reshape(-1, 1)

    if method1 == 3 or method1 == 4:  # save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2
        blocks = [s["param_rs"], s["param_ps"], s["param_as"], s["param_Ks"], s["param_I0s"], s["param_alphas"]]
        labels = ["r", "p", "a", "K0", "I0", "alpha"]
    elif method1 == 5:
        blocks = [s["param_rs"], s["param_ps"], s["param_as"], s["param_Ks"], s["param_I0s"], s["param_alphas"], s["param_ds"]]
        labels = ["r", "p", "a", "K0", "I0", "alpha", "d"]
    else:
        blocks = [s["param_rs"], s["param_ps"], s["param_as"], s["param_Ks"], s["param_I0s"]]
        labels = ["r", "p", "a", "K0", "I0"]

    k = len(labels)

    header = ["time"]
    for label in labels:
        header += _triples(label)
    writetable(names.rolling_csv("parameters-rollingwindow", tstart1, tend1), np.hstack([time] + blocks), header)

    writetable(
        names.rolling_csv("MCSES-rollingwindow", tstart1, tend1),
        np.hstack([time, s["MCSES"][:, :k]]),
        ["time"] + ["%s MCSE" % label for label in labels],
    )

    writetable(
        names.rolling_csv("SCIS-rollingwindow", tstart1, tend1),
        np.hstack([time, s["SCIs"][:, :k]]),
        ["time"] + ["%s SCI" % label for label in labels],
    )


def write_diagnostics_table(names, tstart1, tend1, s):
    """``diagnostics-rollingwindow`` CSV file.

    startAgreement: share of converged optimizer starts that reached the best objective value.
    boot median: median of the bootstrap refits.
    frac at LB / UB: share of bootstrap estimates on the lower / upper bound (NaN = parameter fixed).
    outside CI: 1 when the point estimate lies outside its bootstrap 95% CI.
    best fit on bound: 1 when the point estimate sits on its lower or upper bound.
    """
    labels = ["r", "p", "a", "K0", "I0", "alpha", "d"]
    header = ["time", "startAgreement"]
    for suffix in ("boot median", "frac at LB", "frac at UB", "outside CI", "best fit on bound"):
        header += ["%s %s" % (label, suffix) for label in labels]

    time = np.arange(tstart1, tend1 + 1, dtype=float).reshape(-1, 1)
    table = np.hstack(
        [time, s["startAgreements"], s["param_bootmedians"], s["fracLBs"], s["fracUBs"], s["estimOutsideCIs"], s["bestFitOnBounds"]]
    )
    writetable(names.rolling_csv("diagnostics-rollingwindow", tstart1, tend1), table, header)


def write_aicc_table(names, tstart1, tend1, AICcs):
    """``AICcs-rollingwindow`` CSV file."""
    writetable(
        names.rolling_csv("AICcs-rollingwindow", tstart1, tend1),
        AICcs,
        ["time", "AICc", "AICc part1", "AICc part2", "numparams"],
    )


def write_performance_calibration(names, tstart1, tend1, windowsize1, s):
    """``AICcs-rollingwindow`` and ``performance-calibration`` CSV files."""
    write_aicc_table(names, tstart1, tend1, s["AICcs"])

    time = np.arange(tstart1, tend1 + 1, dtype=float).reshape(-1, 1)
    table = np.hstack([time, np.zeros_like(time) + windowsize1, s["MAECSS"], s["MSECSS"], s["PICSS"], s["WISCSS"]])
    writetable(
        names.rolling_csv("performance-calibration", tstart1, tend1),
        table,
        ["time", "calibration_period", "MAE", "MSE", "Coverage 95%PI", "WIS"],
    )


def write_performance_forecasting(names, tstart1, tend1, forecastingperiod, s):
    """``performance-forecasting`` CSV file."""
    time = np.arange(tstart1, tend1 + 1, dtype=float).reshape(-1, 1)
    table = np.hstack([time, np.zeros_like(time) + forecastingperiod, s["MAEFSS"], s["MSEFSS"], s["PIFSS"], s["WISFSS"]])
    writetable(
        names.rolling_csv("performance-forecasting", tstart1, tend1),
        table,
        ["time", "Horizon", "MAE", "MSE", "Coverage 95%PI", "WIS"],
    )


def write_quantiles(path, quantilesc, quantilesf):
    """The 23 predictive quantiles, calibration rows followed by forecast rows."""
    quantilesc = np.asarray(quantilesc, dtype=float).reshape(-1, len(QUANT_NAMES_RANKED))
    quantilesf = np.asarray(quantilesf, dtype=float).reshape(-1, len(QUANT_NAMES_RANKED))
    writetable(path, np.vstack([quantilesc, quantilesf]), QUANT_NAMES_RANKED)
