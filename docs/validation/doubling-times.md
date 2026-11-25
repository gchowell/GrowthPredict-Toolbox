# Doubling times for windows inside an epidemic: change and test results

Base commit: `b7829ae`, after the other v1.1.0 changes. The comparisons below were run on 25 Nov 2026 with a Python port of `getDoublingTimeCurve.m` (`python/dt_port.py`), not in MATLAB. In MATLAB, the five tests of `tests/test_doublingTime.m` pass, and the doubling times of rows 13–32 before and after the change were reproduced; see [README.md](README.md).

## The problem

`getDoublingTimeCurve.m` follows the cumulative count of the curve it is given and records each time that count doubles. The plotting functions gave it the fitted incidence of the window, so the count restarted at the first observation of every window. For a window taken from the middle of an epidemic the function then reported how fast the window's own running total doubled. For rows 13–32 of the U.S. mpox series, where incidence falls throughout, it reported doublings after 1.1 and 3.1 weeks.

## What changed

| File | Change |
| --- | --- |
| `getDoublingTimeCurve.m` | Optional fourth input `priorCumulative`, the cases counted before the first point of the curve; the count that doubles is `priorCumulative + cumsum(curve)`. Outputs are assigned for a curve of one point. A row vector is accepted. |
| `plotFit_GrowthModels.m`, `plotForecast_GrowthModels.m`, `plotFit_ReproductionNumber.m`, `plotForecast_ReproductionNumber.m` | Pass the sum of the observations before the window. The mean doubling time leaves out bootstrap trajectories without a doubling. |
| `README.md` | Describes the count that doubles. |
| `tests/test_doublingTime.m` | Five tests (new file). |

A window that starts at the first row of the series has no earlier cases, so its results do not change.

## Results (`results/doubling_output.txt`)

### 1. Exponential growth

Incidence `exp(r k)` for k = 0 to 39, windows of 20. The incidence doubles every ln(2)/r.

| r | ln(2)/r | Window start | Count restarted at the window | Count of the series |
| ---: | ---: | ---: | --- | --- |
| 0.1 | 6.93 | 0 | 0.90, 1.59, 2.61, 3.76, 4.86 | same |
| | | 10 | 0.90, 1.59, 2.61, 3.76, 4.86 | 5.10, 5.89, 6.35 |
| | | 20 | 0.90, 1.59, 2.61, 3.76, 4.86 | 6.28, 6.61 |
| 0.2 | 3.47 | 0 | 0.81, 1.34, 1.93, 2.46, 2.88, 3.15, 3.31 | same |
| | | 10 | 0.81, 1.34, 1.93, 2.46, 2.88, 3.15, 3.31 | 3.16, 3.31, 3.40, 3.42, 3.45 |
| | | 20 | 0.81, 1.34, 1.93, 2.46, 2.88, 3.15, 3.31 | 3.40, 3.46, 3.44, 3.47, 3.46 |

- With the restarted count, every window returns the same sequence whatever its position, starting far below ln(2)/r.
- With the count of the series, windows inside the epidemic return values close to ln(2)/r, and they equal the doubling times of the whole series followed from the window start (difference 0).
- The window starting at 0 is unchanged. Its first doublings are short because the series itself starts with its first observation.

### 2. Epidemic that rises and falls

Logistic cumulative curve of 20,000 cases with its midpoint at t = 12, windows of 20.

| Window start | Incidence | Cases so far, and by the window end | Count restarted at the window | Count of the series |
| ---: | --- | --- | --- | --- |
| 6 | rising | 949 → 19,970 | 0.63, 0.96, 1.28, 1.59, 2.42 | 1.43, 1.65, 1.93, 3.31 |
| 12 | falling | 10,000 → 19,999 | 1.00, 2.78 | none |
| 16 | falling | 17,616 → 20,000 | 1.69 | none |
| 20 | falling | 19,640 → 20,000 | 2.01 | none |

### 3. U.S. mpox series, windows of 20

| Rows | Cases before the window | Cumulative count | Curve | Count restarted at the window | Count of the series |
| --- | ---: | --- | --- | --- | --- |
| 1–20 | 0 | 14 → 26,730 | observed | 0.87, 0.89, 0.84, 0.93, 1.06, 1.12, 1.04, 1.17, 1.54, 2.40 | same |
| 4–23 | 64 | 144 → 28,349 | observed | 0.53, 0.81, 0.97, 0.99, 1.04, 1.32, 1.88, 3.81 | 0.96, 1.08, 1.07, 1.04, 1.25, 1.79, 3.22 |
| | | | fitted GLM | 0.45, 0.71, 0.84, 0.91, 1.14, 1.40, 1.99, 3.71 | 0.81, 0.85, 0.98, 1.15, 1.40, 1.88, 3.21 |
| 8–27 | 1,033 | 2,045 → 29,361 | observed | 0.60, 1.02, 1.55, 2.75 | 1.17, 1.67, 2.76 |
| | | | fitted GLM | 0.62, 0.99, 1.52, 2.81 | 1.18, 1.62, 2.82 |
| 13–32 | 11,725 | 14,752 → 29,751 | observed | 1.14, 3.09 | 15.09 |
| | | | fitted GLM | 1.14, 3.11 | 15.55 |

### 4. As the plotting functions report it

Bootstrap of the GLM fit (least squares, normal noise, B = 200); mean and 95% interval of each doubling.

| Rows | Count | Doublings | Mean per trajectory |
| --- | --- | --- | ---: |
| 13–32 | restarted at the window | 1.14 (1.11, 1.17); 3.12 (3.03, 3.22) | 2.13 |
| | of the series | 15.65 (13.65, 17.52) | 15.65 |
| 8–27 | restarted at the window | 0.62 (0.59, 0.65); 0.99 (0.96, 1.02); 1.52 (1.49, 1.55); 2.81 (2.74, 2.88) | 1.48 |
| | of the series | 1.18 (1.13, 1.22); 1.62 (1.58, 1.66); 2.83 (2.76, 2.90) | 1.87 |

Every trajectory reached each of these doublings. For rows 13–32 the values with the restarted count agree with the MATLAB output from before the change: 1.136 (1.11, 1.17), 3.112 (3.03, 3.20) and a mean of 2.124. With the count of the series, MATLAB (B = 300) reports one doubling after 15.61 weeks (13.44, 17.32), reached by every trajectory.

## Limits

- **Start of the series.** The first doublings of a series are short when it begins after the outbreak did, because the earlier cases are unknown.
- **Series given as cumulative counts** are converted to incidence by the runners and plotting functions before this step, so the count before the window is still the sum of the earlier rows.
- **No doubling.** A window whose count never doubles gives one row of `NaN` with fraction 0 in `doublingtimes-*.csv`, and a mean doubling time of `NaN`.

## Reproducing

```bash
pip install numpy scipy numba
python docs/validation/python/test_doubling.py /tmp/doubling.json "forecasting_growthmodels code/input"
```

In MATLAB, from the code folder: `results = runtests('tests/test_doublingTime.m')`.
