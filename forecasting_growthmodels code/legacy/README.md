# Legacy scripts

These files come from earlier analyses. The runners, plotting functions and tests do not call them, and they are not on the MATLAB path when you work from the code folder. Most read data or result files that are not in the repository, so they are kept for reference only.

| File | What it was |
|---|---|
| `temp_Rt_slidingwindow.m`, `temp_Rt_slidingwindow2.m`, `temp_Rt_slidingwindowE5.m` | Sliding-window reproduction-number experiments on simulated SEIR epidemics, compared with EpiEstim estimates. They read CSV files that are not in the repository. |
| `temp_Rt_slidingwindow_cholera.m`, `temp_Rt_slidingwindow_choleraDistrict.m` | The same sliding-window analysis for the cholera series and for cholera district epidemics. |
| `plotRt_cholera.m` | Collects the last row of `Rt-*.csv` from 80 cholera windows in `output/` and plots them. |
| `Read_Rs_choleradistricts.m` | Gathers saved reproduction numbers for 28 cholera districts from earlier MAT files in `output/`. |
| `estimateGrowthRatesCholera.m` | Loops over the rows of `cholera-data-Sally.txt` (not in the repository) and fits the growth phase of each. **It overwrites `input/cholera-epidemic.txt`.** |
| `plotParamsByErrorStructure_simulated1.m` | Compares parameter estimates across error structures, read from earlier output files for simulated data. |
| `timingStudy.m` | Times one call of `Run_Forecasting_GrowthModels`. |
| `get_Rt_old.m` | `get_Rt.m` as it was before v1.1.0, with the errors corrected in that version. It declares a function named `get_Rt`. |
| `options_Rt_xinyi.m` | Generation-interval settings from an earlier analysis (gamma, mean 4, variance 8). It declares a function named `options_Rt`. |

To run one from an old script, call `addpath('legacy')` from the code folder first. Do not add the folder to the path otherwise: `get_Rt_old.m` and `options_Rt_xinyi.m` declare the names of current functions.
