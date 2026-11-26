# Generation-interval defaults: change and test results

Base commit: `b7829ae`, after the `r`-interval and `get_Rt.m` changes. The comparisons below were run on 25 Nov 2026 with the Python port (`python/test_gi_defaults.py`), not in MATLAB. In MATLAB, the six tests of `tests/test_optionsRt.m` pass and the value at the end of rows 13–32 was reproduced; see [README.md](README.md).

## What changed

| File | Before | After |
| --- | --- | --- |
| `options_Rt.m` | Mean 5 days, SD 8 days (cholera) | Mean 8.5 days, SD 5.0 days (mpox, U.S. case pairs, symptom onset), with the source in the header |
| `options_Rt_mpox.m` | Mean 3.8, variance 2·(1/3.8)² | Same values as `options_Rt.m` |
| `options_Rt_cholera.m` | Mean 5 days, SD 8 days; function named `options_Rt` | Same values; function name matches the file; source and the missing one-day shift noted |
| `README.md` | No statement of the default | One sentence naming the default, its source and the earlier default |
| `tests/test_optionsRt.m` | | Six tests (new file) |
| `tests/test_getRt.m` | Nine tests | Eight: the test that pinned the old `options_Rt_mpox.m` values moved to the new file |

All values are in weeks, because the default series is weekly. `get_Rt.m` is unchanged.

Sources: Madewell ZJ, et al. Emerg Infect Dis. 2023;29(4), doi:10.3201/eid2904.221622 (mpox). The cholera values follow the Yemen 2016–18 analysis, which used a gamma distribution shifted by one day with mean 5 days and SD 8 days. That analysis is Camacho A, et al. Lancet Glob Health. 2018;6(6):e680-e690, doi:10.1016/S2214-109X(18)30230-4; the interval is given in its appendix.

## The distributions (weekly units)

| | Mean | Variance | Shape | Mass within 1 day | 95% of mass within |
| --- | ---: | ---: | ---: | ---: | ---: |
| Former default (cholera) | 0.714 | 1.306 | 0.39 | 40.7% | 3.0 weeks |
| New default (mpox, symptom onset) | 1.214 | 0.510 | 2.89 | 0.7% | 2.6 weeks |
| Mpox, rash onset (7.0 d, SD 4.2 d) | 1.000 | 0.360 | 2.78 | 1.3% | 2.1 weeks |

## Effect on the bundled U.S. mpox example

Rt of the fitted GLM curve (least squares, `fixI0 = 1`).

| Rows | Generation interval | 5 weeks in | 10 weeks in | Window end |
| --- | --- | ---: | ---: | ---: |
| 1–20 | Former default | 1.31 | 1.15 | 0.754 |
| | New default | 1.84 | 1.30 | 0.654 |
| | Rash onset | 1.65 | 1.23 | 0.702 |
| 8–27 | Former default | 1.04 | 0.86 | 0.755 |
| | New default | 1.03 | 0.77 | 0.687 |
| | Rash onset | 1.01 | 0.80 | 0.735 |
| 13–32 | Former default | 0.88 | 0.81 | 0.784 |
| | New default | 0.77 | 0.73 | 0.720 |
| | Rash onset | 0.81 | 0.77 | 0.764 |

In rows 1–20 the time at which Rt falls below 1 moves from t = 13.0 to t = 12.7.

## Check against the exact value

Exponential growth at r = 0.3 per week with the new default: exact R = 1.4093, estimate at the window end 1.4004 (−0.63%). At r = −0.2: exact 0.7759, estimate 0.7725 (−0.44%).

## Notes

- `docs/images/Rt.png` depends on the default; regenerate it with `make_readme_figures` after changing `options_Rt.m`.
- Rt results produced with the default options file change; see `CHANGELOG.md`.
- `options_Rt_xinyi.m` (mean 4, variance 8) is untouched. (In v1.1.1 it moved to the code folder's `legacy` folder.)

## Reproducing

```bash
python docs/validation/python/test_gi_defaults.py /tmp/gi_defaults.json "forecasting_growthmodels code/input"
```

In MATLAB, from the code folder: `results = runtests('tests/test_optionsRt.m')`.
