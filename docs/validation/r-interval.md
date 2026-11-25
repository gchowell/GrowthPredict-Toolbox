# Interval of `r` in `fit_model.m`: change and test results

Base commit: `b7829ae` (22 Nov 2026). The comparisons below were run on 24 Nov 2026 with a Python port of the fit (`python/gp_port.py`), not in MATLAB. In MATLAB, the seven tests of `tests/test_rInterval.m` pass and the fits of rows 1–20 and 13–32 were reproduced; see [README.md](README.md).

## What changed

| File | Change |
| --- | --- |
| `fit_model.m` | GGM, GLM and GRM: `r` in `[1e-4, 10*max(y)]/DT` instead of the slope-based interval. New option `fitopts.rBounds = [lower upper]`. Header documents both. |
| `Run_Fit_GrowthModels.m`, `Run_Forecasting_GrowthModels.m` | Optional last argument: a struct merged into the optimizer settings, e.g. `struct('rBounds',[1e-3 50])`. |
| `README.md` | One paragraph in "Optimization and bootstrap" describes the bounds of `r`. |
| `tests/test_rInterval.m` | Seven regression tests (new file). |

EXP, LM, RICH and GOM are unchanged: with `p` fixed at 1, `r` is a rate per unit time and the slope-based interval has the right units.

## What the Python port is

`python/gp_port.py` follows `fit_model.m` for the models, the incidence curve, the least-squares and Poisson objectives, the parameter box, the seed, the start points (10 Latin-hypercube points, 3 jittered seeds, the seed), the unit-box and log10 coordinates, and the retry with a wider `r` interval. It uses `scipy.optimize.least_squares` in place of `fmincon`/`MultiStart`, so optimizer behaviour (start agreement, timing) is indicative only.

Checks of the port: its ODE solver agrees with LSODA to 7e-10 (relative) over 1,000 random parameter sets, and it reproduces the README quick start (r = 1.73, p = 0.853, K = 28,547).

## 1. All bundled series (`results/matrix_summary.txt`)

18 series, 149 windows (rolling windows of 20, or 12 for the 16-point HIV series, plus each full series), GGM/GLM/GRM, least squares and Poisson, `fixI0 = 1`: 894 fits. "Best known" is the lowest objective over five runs per case (old box, new box, an initial-rate parametrization, and 40-start runs of the last two).

| | Old interval | New interval |
| --- | ---: | ---: |
| `r` on a bound | 450 (50.3%) | 28 (3.1%), all GRM |
| Reaches best-known objective | 51.3% | 99.1% |
| `p` on a bound | 18.7% | 44.4% |
| Any of r, p, a, K on a bound | 78.7% | 58.8% |
| Start agreement below 0.5 | 2.5% | 2.7% |
| Seconds per fit (port) | 0.56 | 0.81 |

- New versus old objective: lower in 435 fits, equal in 458, higher in 1 (GRM, +0.2%, an optimizer miss).
- Of the 450 fits with `r` on the old bound (always the upper one), the new interval improved 411.
- Least squares, where improved (226 fits): SSE ratio old/new has median 1.67; above 2 in 40%, above 10 in 12%.
- `p` differed by more than 0.1 in 87% of the improved fits (median difference 0.35).
- GGM and GLM: the largest fitted `r` was 2.18 times `max(y)`, against a bound of 10 times.
- Four-step point forecasts (654 fits with held-out data): changed in 328; error lower in 229, higher in 99; median ratio new/old 0.85.
- `fixI0 = 0` (GLM, least squares, 66 fits): new interval reached the best-known objective in all; old in half.

## 2. Simulated Poisson data, 60 replicates each (`results/recovery.json`)

| Scenario | Interval | Median `r` | Bias of `p` | `r` on bound | 7-step error / mean |
| --- | --- | ---: | ---: | ---: | ---: |
| A: GLM r=0.18, p=0.94, K=1e4 | old | 0.181 | 0.000 | 0% | 0.076 |
| | new | 0.181 | 0.000 | 0% | 0.076 |
| B: GLM r=40, p=0.40, K=5e5 | old | 10 | +0.171 | 100% | 0.258 |
| | new | 40.4 | 0.000 | 0% | 0.012 |
| C: GGM r=25, p=0.30 | old | 10 | +0.124 | 100% | 0.209 |
| | new | 24.8 | +0.001 | 0% | 0.018 |
| D: GLM r=3, p=0.70, K=2e4 | old | 2.99 | 0.000 | 0% | 0.028 |
| | new | 2.99 | 0.000 | 0% | 0.028 |

In B the old interval also returned a median K of 1.0e5 against a true 5e5; the new one returned 4.9e5.

## 3. Bootstrap on the mpox series, GLM, least squares, B = 200 (`results/bootstrap.json`)

| Rows | Interval | `r` (95% interval) | Refits with `r` at cap | `p` (95% interval) | 4-step MAE |
| --- | --- | --- | ---: | --- | ---: |
| 1–20 | old and new | 1.73 (1.49, 2.08) | 0% | 0.853 (0.828, 0.874) | 193 |
| 8–27 | old | 10 (6.05, 10) | 58% | 0.672 (0.667, 0.725) | 52 |
| 8–27 | new | 27.0 (21.7, 35.2) | 0% | 0.558 (0.527, 0.583) | 27 |
| 13–32 | old | 10 (2.13, 10) | 53% | 0.676 (0.664, 0.856) | no held-out data |
| 13–32 | new | 803 (545, 1238) | 0% | 0.177 (0.128, 0.223) | no held-out data |

## What the change does not fix

- **GRM ridge.** In 28 of 298 GRM fits `r` still reached `10*max(y)`, all with `a` below 0.2, because `r` and `a` trade off as `a` approaches 0. Six of them were above the best-known objective (largest gap 25%, the rest at most 5.5%). A variant with `a >= 0.05` and a 20-fold wider `r` interval moved the flag from `r` to `a` and made 4 of 40 other fits marginally worse, so it was not adopted (`results/grm_variant_summary.txt`).
- **`p` at 0.01.** With `r` free, `p` reached its lower bound in 30% of fits (6% before), 59% for GGM. These are windows that level off or decline, which the model cannot follow; the old cap on `r` hid this.
- **K at its lower bound.** In the bootstrap above, K sat on its lower bound (the window's cumulative count) in 28–36% of refits.
- **Initial-rate parametrization.** Searching on `i0 = r*C0^p*(1 - C0/K)` with the box `[yref/20, 20*yref]` was faster (0.39 s per fit) but that box bound in 74 fits and missed the best-known objective in 75. It needs a wider box before it can replace this change.

## Data notes

- `cholera-district-epidemic.txt` stacks several districts (the time index resets), so its windows were excluded. `fit_model.m` does not check that the time column increases.
- 65 windows with fewer than 10 events or fewer than 4 non-zero counts were skipped.

## Reproducing

From the repository root:

```bash
pip install numpy scipy numba
IN="forecasting_growthmodels code/input"; V=docs/validation
python $V/python/validate_port.py "$IN/Most_Recent_Timeseries_US-CDC.txt"
python $V/python/run_matrix.py "$IN" $V/results/matrix.jsonl        # about 40 minutes on 2 cores
python $V/python/analyze_matrix.py $V/results/matrix.jsonl
python $V/python/run_recovery.py $V/results/recovery.json 60
python $V/python/run_bootstrap.py "$IN/Most_Recent_Timeseries_US-CDC.txt" 200 $V/results/bootstrap.json
python $V/python/run_grm_variant.py "$IN" $V/results
```

`matrix.jsonl` (one line per fit, 3 MB) is not kept in the repository; `run_matrix.py` regenerates it and `results/matrix_summary.txt` is its summary.

In MATLAB, from the code folder: `results = runtests('tests/test_rInterval.m')`.
