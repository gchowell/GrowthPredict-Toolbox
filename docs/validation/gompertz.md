# Gompertz model, time origin: change and test results

Base commit: `b7829ae`, after the other v1.1.0 changes. The comparisons below were run on 25 Nov 2026 with the Python port of the fit (`python/gp_port.py`), not in MATLAB. In MATLAB, the six tests of `tests/test_gompertz.m` pass and the fit of rows 13–32 was reproduced (r = 0.6732, a = 0.3939, sum of squares 4.525e5); see [README.md](README.md).

## The problem

The Gompertz model `dC/dt = r*C*exp(-a*t)` is the only model in the toolbox whose growth rate depends on time. The solver used the time labels of the window as `t`. The same 20 observations labelled 12 to 31 instead of 0 to 19 therefore needed an `r` larger by `exp(12*a)` to give the same curve. Two things followed:

- **`r` had no fixed meaning.** It was the growth rate extrapolated back to time label 0, not the growth rate at the start of the window.
- **The fit was limited for later windows.** The interval of `r` is about 0.1 to 10, so once `exp(a*t0)` exceeded it, `r` stopped at its upper bound and the curve could not follow the data.

## What changed

| File | Change |
| --- | --- |
| `growthIncidenceCurve.m` | For the Gompertz model, time is measured from the first point of the curve. The other models do not depend on `t` and are solved as before. |
| `fit_model.m` | The Gompertz bounds and seed, which `initialParams.m` gives per observation step, are divided by `DT`. No effect with `DT = 1`, the value the runners use. |
| `modifiedLogisticGrowth.m`, `README.md` | State that `t` is the time since the first observation of the window. |
| `tests/test_gompertz.m` | Six tests (new file). |

A window whose first time label is 0 gives the same result as before.

## What was tested and not adopted

The review that led to this change also said that the bounds of `r` ignore the time unit. The tests do not support changing them:

- **`r`.** With time measured from the window start, `r` never reached a bound in 298 fits; the estimates ran from 0.14 to 7.9 against an interval of about 0.1 to 10. Because the cumulative count starts at the first observation, `r` is a rate per observation step and falls in this range for daily, weekly and yearly series alike.
- **`a`.** A wider interval (`r` and `a` in `[1e-4, 10]`) lowered the objective in 17 of 298 fits, by a median of 16% in sum of squares, by letting `a` go to zero, which turns the model into exponential growth. The four-step forecast became worse in 15 of those 17 (median error 58% higher). The current lower bound of `a`, one tenth of the guess of `initialParams.m`, keeps the implied final size finite and was left as it is.

## Results

### 1. All bundled series (`results/gompertz_summary.txt`)

149 windows of 18 series, least squares and Poisson, `fixI0 = 1`: 298 Gompertz fits.

| | Before | After |
| --- | ---: | ---: |
| `r` on its upper bound | 153 (51.3%) | 0 |
| `a` on its lower bound | 19 | 17 |
| Reaches best-known objective | 40.9% | 94.0% |
| Start agreement below 0.5 | 34.2% | 0.3% |

- **Objective:** lower in 159 fits, equal in 139, higher in none. In the least-squares fits that improved, the sum of squares fell by a median factor of 3.65, more than 2 in 74% and more than 10 in 16%.
- **Where:** all 159 are windows whose first time label is not 0 (240 fits). The 58 fits whose first label is 0 are unchanged.
- **Four-step point forecasts** (218 fits with held-out data): changed in 118; error lower in 96, higher in 22; median ratio after/before 0.38.

### 2. Same observations, different time labels (`results/gompertz_tests_output.txt`)

U.S. mpox rows 1–20, least squares.

| First label | Before: `r`, `a`, sum of squares | After: `r`, `a`, sum of squares |
| ---: | --- | --- |
| 0 | 1.474, 0.188, 2.72e6 | 1.474, 0.188, 2.72e6 |
| 1 | 1.778, 0.188, 2.72e6 | same as label 0 |
| 7 | 5.481, 0.188, 2.72e6 | same as label 0 |
| 12 | 10 (bound), 0.168, 3.94e6 | same as label 0 |
| 30 | 10 (bound), 0.086, 2.66e7 | same as label 0 |
| 100 | 10 (bound), 0.031, 4.11e7 | same as label 0 |

### 3. Time step

Same observations with `DT` of 1, 7, 1/7 and 30. After the change the curve is identical and `r*DT` and `a*DT` are 1.474 and 0.188 in each case. Before it, `DT = 1/7` ended on a bound and `DT = 30` gave a sum of squares of 5.8e7 against 2.7e6.

### 4. Simulated Gompertz epidemics

True `r = 0.8`, `a = 0.15`, 25 points with Poisson noise, 60 replicates, Poisson maximum likelihood, seven-step forecast.

| First label | | Reported `r` (median) | Growth rate at the first observation | `a` (median) | `r` on bound | Forecast error / mean |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | before and after | 0.798 | 0.798 | 0.150 | 0% | 0.040 |
| 10 | before | 3.569 | 0.798 | 0.150 | 0% | 0.040 |
| | after | 0.798 | 0.798 | 0.150 | 0% | 0.040 |
| 40 | before | 9.954 | 0.460 | 0.077 | 100% | 3.191 |
| | after | 0.798 | 0.798 | 0.150 | 0% | 0.040 |

With a first label of 10 the curve was right before the change but `r` was reported as 3.57 for a true growth rate of 0.8.

### 5. Bootstrap on the mpox series, least squares, B = 200

| Rows | | RMSE | `r` (95% interval) | Refits with `r` at its upper bound | `a` (95% interval) |
| --- | --- | ---: | --- | ---: | --- |
| 8–27 | before | 96.4 | 8.22 (7.58, 8.84) | 0% | 0.301 (0.293, 0.308) |
| | after | 96.4 | 0.998 (0.974, 1.02) | 0% | 0.301 (0.293, 0.308) |
| 13–32 | before | 374.8 | 8.33 (3.42, 8.33) | 54% | 0.245 (0.190, 0.251) |
| | after | 150.4 | 0.673 (0.625, 0.720) | 0% | 0.394 (0.359, 0.427) |

## Limits

- **`a` on its lower bound.** In 17 of 298 fits, windows of early growth with no visible slowing, `a` stays on its lower bound. The data do not determine when growth slows there; the runners report the bound.
- **`K` is not a Gompertz parameter.** Parameter exports show `K = 1` for this model. The final size implied by a fit is `C0*exp(r/a)`.
- **Other models.** Rolling windows of the other models reset the cumulative count at each window but do not depend on the time labels.

## Reproducing

```bash
pip install numpy scipy numba
IN="forecasting_growthmodels code/input"; V=docs/validation
python $V/python/run_gompertz.py "$IN" /tmp/gompertz.jsonl        # about 2 minutes on 2 cores
python $V/python/analyze_gompertz.py /tmp/gompertz.jsonl
python $V/python/test_gompertz.py /tmp/gompertz_tests.json "$IN"
```

In MATLAB, from the code folder: `results = runtests('tests/test_gompertz.m')`.
