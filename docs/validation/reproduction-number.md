# Reproduction number (`get_Rt.m`): change and test results

Base commit: `b7829ae`. The comparisons below were run on 25 Nov 2026 with Python ports of the old and new function (`python/rt_port.py`), not in MATLAB. In MATLAB, the tests of `tests/test_getRt.m` pass and the direct checks were reproduced; see [README.md](README.md).

## What changed

| File | Change |
| --- | --- |
| `get_Rt.m` | Generation interval measured as elapsed time from the start of the curve. Fixed interval placed at exactly `g` (rounded to the nearest 0.01) with a rounded index. Input checks. Optional third output `giCoverage`. |
| `plotFit_ReproductionNumber.m`, `plotForecast_ReproductionNumber.m` | `Rt-*.csv` rows are labelled with the time each estimate refers to: one row per unit of time after the window start, ending at the last time. |
| `options_Rt_mpox.m` | `1//3.8` corrected to `1/3.8`; function name matches the file. |
| `README.md` | Rt paragraphs describe the new behaviour and say which earlier results to rerun. |
| `tests/test_getRt.m` | Regression tests (new file): nine at first, eight after the option-file test moved to `tests/test_optionsRt.m`. |

## Results (`results/rt_tests_output.txt`)

### 1. Exponential growth, where R is known exactly

Weekly samples, window of 20, estimate at the window end.

| Generation interval | r | Exact R | Old, start 0 | Old, start 1 | Old, start 12 | New, any start |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Gamma, mean 5/7, variance (8/7)² (default) | 0.3 | 1.1863 | 1.1822 | 6.06 | 8,227 | 1.1822 |
| | 0.1 | 1.0678 | 1.0677 | 4.97 | 6,412 | 1.0677 |
| | −0.2 | 0.8371 | 0.8341 | 3.17 | 3,596 | 0.8341 |
| Gamma, mean 1.4, SD 0.6 | 0.3 | 1.4987 | 1.4899 | 1.65 | 3.4e14 | 1.4899 |
| Exponential, mean 1 | 0.3 | 1.3000 | 1.2926 | 3.51 | 2.1e5 | 1.2926 |
| Fixed, 1 | 0.3 | 1.3499 | 1.3534 | 1.3534 | 1.3534 | 1.3499 |

- New function: largest relative error over all 36 combinations is 0.59%, which comes from interpolating weekly samples linearly.
- Old function at start 1 or 12: median relative error 225%.
- For a curve starting at zero, the gamma and exponential results are identical before and after. Only the fixed interval changes there (by one grid step).

### 2. Moving the window start

The same curve labelled from 0, 1, 5, 12, 50 or 1,000 gives the same values with the new function (largest relative change 2.5e-14). With the old function the gamma and exponential values changed by factors of 1e5 to 1e14.

### 3. Fixed interval over every lag 0.01 to 10.00

| | Old | New |
| --- | --- | --- |
| Lags that stop with an indexing error | 158 of 1,000 (0.14, 0.28, 0.29, 0.47, 0.56, 0.57, ...) | 0 |
| Lag actually applied | g + 0.01 | g (error below 1e-14) |

A lag off the grid is rounded: 5/7 is applied as 0.71. A lag longer than the curve gives all `NaN`; a lag below 0.005 is rejected.

The count of 158 uses lags typed as two-decimal numbers. The 148 quoted in the review came from computing the lags as k × 0.01.

### 4. Renewal simulation with a changing R(t)

Incidence simulated from the renewal equation with R(t) falling from 2.6 to 0.7 (gamma, mean 1.4, SD 0.6), then sampled at integer times.

| Window | True R at end | Old | New | New: largest error after 6 units of history |
| --- | ---: | ---: | ---: | ---: |
| 0–19 | 0.881 | 0.886 | 0.886 | 4.1% |
| 5–24 | 0.721 | 18,020 | 0.719 | 2.8% |
| 12–31 | 0.701 | 2.9e14 | 0.697 | 0.7% |
| 20–39 | 0.700 | NaN | 0.697 | 0.5% |

### 5. Bundled series, Rt of the fitted GLM curve at the window end (default generation interval)

| Series | Rows | First time label | Old | New |
| --- | --- | ---: | ---: | ---: |
| `Most_Recent_Timeseries_US-CDC` | 1–20 | 0 | 0.754 | 0.754 |
| | 4–23 | 3 | 10.9 | 0.730 |
| | 8–27 | 7 | 146 | 0.755 |
| | 13–32 | 12 | 3,131 | 0.784 |
| `weekly-Mpox-cases-Canada-05-31-2023` | 1–20 | 1 | 3.35 | 0.857 |
| `weekly-Mpox-cases-US_CDC-05-24-2023` | 1–20 | 1 | 3.04 | 0.804 |

The third series with a first time label of 1 (the stochastic SIR simulation) has a fitted curve that is flat at zero by the window end, so its Rt is not informative under either version.

### 6. CSV export

For a window labelled 0 to 19, the old file had 20 rows labelled 0.01, 1.01, ..., 18.01, 19 with estimates for 0.02, 1.02, ..., 18.02, 19. The new file has 19 rows labelled 1, 2, ..., 19 with the estimates for those times. A later change leaves out the first rows; see [reproduction-number-coverage.md](reproduction-number-coverage.md). `plotRt.m` reads the last row, which is the window end in both.

## What the change does not fix

- **Start of every window.** Incidence before the window is not used, so early estimates are too high. For exponential growth at r = 0.2 with the generation interval that was then the default (mean 5 days, SD 8 days) the error is +20% one time unit in, +6% at two and under 1% from about four. A later change uses `giCoverage` to leave these estimates out of the plots and the CSV; see [reproduction-number-coverage.md](reproduction-number-coverage.md).
- **Generation-interval values.** This change left the default (mean 5 days, SD 8 days) and the values in `options_Rt_mpox.m` as they were, apart from the syntax error. Both were changed afterwards; see [generation-interval-defaults.md](generation-interval-defaults.md).
- **`get_Rt_old.m`** is still in the folder and still declares a function named `get_Rt`. (In v1.1.1 it moved to the code folder's `legacy` folder.)

## Reproducing

```bash
pip install numpy scipy numba
python docs/validation/python/test_rt.py /tmp/rt_tests.json "forecasting_growthmodels code/input"
```

In MATLAB, from the code folder: `results = runtests('tests/test_getRt.m')`.
