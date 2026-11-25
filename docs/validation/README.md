# Validation notes for v1.1.0

Each change in v1.1.0 was checked against an independent Python re-implementation of the relevant MATLAB function, and then confirmed in MATLAB through regression tests and direct runs. The notes give the comparisons; `python/` holds the re-implementations and scripts; `results/` holds their output.

| Change | Note | MATLAB tests |
| --- | --- | --- |
| Interval of `r` for GGM, GLM and GRM | [r-interval.md](r-interval.md) | `tests/test_rInterval.m` (7) |
| Renewal calculation in `get_Rt.m` | [reproduction-number.md](reproduction-number.md) | `tests/test_getRt.m` (8) |
| Default generation interval | [generation-interval-defaults.md](generation-interval-defaults.md) | `tests/test_optionsRt.m` (6) |
| Reproduction number at the start of a window | [reproduction-number-coverage.md](reproduction-number-coverage.md) | `tests/test_selectRtRows.m` (7) |
| Doubling times for windows inside an epidemic | [doubling-times.md](doubling-times.md) | `tests/test_doublingTime.m` (5) |
| Time origin of the Gompertz model | [gompertz.md](gompertz.md) | `tests/test_gompertz.m` (6) |
| Titles of the parameter histograms | [figure-titles.md](figure-titles.md) | `tests/test_parameterTitle.m` (10) |

## What was run in MATLAB

Run by the maintainer on 25 Nov 2026 from `forecasting_growthmodels code`, on the U.S. mpox series (`Most_Recent_Timeseries_US-CDC.txt`).

- **Regression tests.** `runtests('tests')` passed 49 of 49 across the seven files in the table.
- **Direct checks.** The Python value is from the re-implementations in `python/`.

| Check | Python | MATLAB |
| --- | --- | --- |
| GLM fit, rows 13–32: r, p, K, sum of squares | 803.1, 0.177, 18,153, 2.3225e4 | 801.8, 0.178, 18,153, 2.323e4 |
| Same fit with `rBounds = [1e-3 10]` | r = 10, p = 0.676, 9.51e5 | r = 10, p = 0.676, 9.51e5 |
| Quick start (rows 1–20), bootstrap: r | 1.73 (1.49, 2.08) | 1.7 (1.5, 2.1) |
| Quick start, bootstrap: p | 0.853 (0.828, 0.874) | 0.85 (0.83, 0.87) |
| Rows 13–32, bootstrap: r | 803 (545, 1238) | 8e2 (5.7e2, 1.2e3) |
| Rt, exponential growth r = 0.3, cholera interval, window starting at 0 and at 12 | 1.1822 and 1.1822 | 1.1822 and 1.1822 |
| Rt, same curve, fixed interval of 1 | 1.349859 | 1.349859 |
| Rt of rows 13–32 at t = 13, before the coverage filter | 2.04 | 2.045 |
| Rt of rows 13–32, first exported row (t = 15) | 0.85, coverage 0.977 | 0.850 (0.837, 0.863), coverage 0.9767 |
| Rt of rows 13–32, last row (t = 31) | 0.720 | 0.721 (0.706, 0.735) |
| Rows exported for rows 13–32 | 17 | 17 |
| Doubling times of rows 13–32, count restarted at the window | 1.14 (1.11, 1.17); 3.12 (3.03, 3.22) | 1.136 (1.11, 1.17); 3.112 (3.03, 3.20) |
| Doubling time of rows 13–32, count of the series | 15.65 (13.65, 17.52), reached by 100% | 15.61 (13.44, 17.32), reached by 100% |
| Gompertz fit, rows 13–32: r, a, sum of squares | 0.6732, 0.3939, 4.525e5 | 0.6732, 0.3939, 4.525e5 |

The Python fit uses the fitted curve where MATLAB reports the median of bootstrap curves, and SciPy in place of `fmincon`, so small differences are expected.

## Running the Python scripts

From the repository root, with Python 3 and `pip install numpy scipy numba`. Each note ends with its own commands. `python/gp_port.py` is the port of `fit_model.m`, `python/rt_port.py` of `get_Rt.m` and `selectRtRows.m`, and `python/dt_port.py` of `getDoublingTimeCurve.m`.

The file with one line per fit of the comparison in `r-interval.md` (`matrix.jsonl`, 3 MB) is not kept here; `python/run_matrix.py` regenerates it in about 40 minutes on two cores.
