# Rt at the start of a window: change and test results

Base commit: `b7829ae`, after the `r`-interval, `get_Rt.m` and generation-interval changes. The comparisons below were run on 25 Nov 2026 with the Python ports (`python/test_coverage.py`), not in MATLAB. In MATLAB, the seven tests of `tests/test_selectRtRows.m` pass and the exported rows for rows 13–32 were reproduced; see [README.md](README.md).

## What changed

| File | Change |
| --- | --- |
| `selectRtRows.m` (new) | From the coverage returned by `get_Rt.m`, picks the estimates whose incidence history covers at least 95% of the generation-interval distribution. |
| `plotFit_ReproductionNumber.m`, `plotForecast_ReproductionNumber.m` | The Rt curve starts where coverage reaches 95%. `Rt-*.csv` omits the rows below it and gains a `GI coverage` column. The y-axis is scaled to the curve shown (it was the mid-window value plus 5). A stray console print is removed, and the Rt = 1 line starts at the window start. |
| `README.md` | Describes the filter, the new column and its limit. |
| `tests/test_selectRtRows.m` | Seven tests (new file). |

`plotRt.m` and `plotRt_cholera.m` read the last row and the first four columns of each `Rt-*.csv`, so they are unaffected. (In v1.1.1 `plotRt_cholera.m` moved to the code folder's `legacy` folder.) The threshold is the constant `minCoverageRt = 0.95` in each plotter.

## Results (`results/coverage_output.txt`)

### 1. Rows kept for a window labelled 12 to 31

| Generation interval | Reliable from | CSV rows | Dropped |
| --- | ---: | --- | --- |
| Mpox default (8.5 d, SD 5.0 d) | t = 14.58 | 15 to 31 (17 rows) | 13, 14 |
| Mpox, rash onset (7.0 d, SD 4.2 d) | t = 14.15 | 15 to 31 | 13, 14 |
| Cholera (5 d, SD 8 d) | t = 15.00 | 15 to 31 | 13, 14 |
| Exponential, mean 1 | t = 15.00 | 15 to 31 | 13, 14 |
| Fixed, 1 | t = 13.00 | 13 to 31 (19 rows) | none |

### 2. Error of the dropped and kept rows where R is known exactly

Exponential growth or decline at rate r per week; error relative to the exact R.

| Generation interval | r | Dropped rows | Largest error among kept rows |
| --- | ---: | --- | ---: |
| Mpox default | −0.3 | +167%, +23% | 3.9% |
| | −0.2 | +150%, +20% | 3.4% |
| | +0.1 | +109%, +13% | 1.8% |
| | +0.3 | +88%, +9% | 0.6% |
| Mpox, rash onset | −0.3 | +96%, +11% | 0.9% |
| Cholera | −0.3 | +63%, +32% | 19.4% |
| | −0.2 | +47%, +22% | 11.9% |
| | +0.3 | +16%, +5% | 1.3% |
| Exponential, mean 1 | −0.3 | +97%, +32% | 12.9% |
| | +0.3 | +37%, +7% | 1.5% |

### 3. Effect of the threshold

Largest error among kept rows, and the first exported row.

| Generation interval | Threshold | First row | r = −0.3 | r = −0.2 | r = +0.1 | r = +0.3 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Mpox default | 0.95 | 15 | 3.9% | 3.4% | 1.8% | 0.6% |
| | 0.99 | 16 | 0.9% | 0.4% | 0.2% | 0.6% |
| Cholera | 0.95 | 15 | 19.4% | 11.9% | 3.5% | 1.3% |
| | 0.99 | 18 | 5.7% | 2.6% | 0.3% | 0.3% |
| Exponential, mean 1 | 0.95 | 15 | 12.9% | 9.5% | 3.8% | 1.5% |
| | 0.99 | 17 | 2.2% | 1.4% | 0.4% | 0.6% |

A threshold of 0.90 gives the same rows as 0.95 for these three intervals.

### 4. U.S. mpox series, Rt of the fitted GLM curve, mpox default

| Rows | Dropped | First kept | Window end | Largest kept Rt |
| --- | --- | --- | ---: | ---: |
| 1–20 | t = 1: 3.20; t = 2: 2.33 | t = 3: 2.12 | 0.654 | 2.12 |
| 8–27 | t = 8: 2.91; t = 9: 1.63 | t = 10: 1.31 | 0.687 | 1.31 |
| 13–32 | t = 13: 2.04; t = 14: 1.02 | t = 15: 0.85 | 0.720 | 0.85 |

The value at t = 13 for rows 13–32 (2.04) matches the 2.045 in the MATLAB export from before this change.

### 5. Window too short

A generation interval with mean 10 and SD 5 on a window of 5 time units reaches a coverage of 0.14. No row is reliable; the CSV holds the last time alone with its coverage, the plotters warn, and the whole curve is drawn.

## Limits

- **Heavy-tailed intervals with falling incidence.** At 95% coverage the first kept row can still be 12% to 19% too high for the cholera and exponential intervals. The coverage column allows a stricter cut; 0.99 brings these to 6% or less.
- **Rows 1–20 start at the beginning of the series,** so there is little earlier incidence to miss; the filter still drops t = 1 and t = 2 because it cannot tell a series start from a window cut out of a longer series.
- **The threshold is a constant** in the two plotters, not an option.

## Reproducing

```bash
python docs/validation/python/test_coverage.py /tmp/coverage.json "forecasting_growthmodels code/input"
```

In MATLAB, from the code folder: `results = runtests('tests/test_selectRtRows.m')`.
