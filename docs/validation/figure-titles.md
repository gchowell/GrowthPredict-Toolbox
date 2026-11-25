# Titles of the parameter histograms: change and test results

Base commit: `664ff0f`, after the other v1.1.0 changes. The checks below were run on 25 Nov 2026 with a Python port of the title text and font-size rule (`python/title_port.py`) and text rendered by matplotlib, not in MATLAB. In MATLAB, the ten tests of `tests/test_parameterTitle.m` pass, including those that measure the rendered titles, and `docs/images/fit.png` and `forecast.png` were regenerated with the new titles; see [README.md](README.md).

## The problem

The top row of the fit and forecast figures has four histograms (`r`, `p`, `a`, `K`), drawn with `subplot(2,4,k)`. Each had a one-line title such as `K=2.9e+04 (95% CI:2.7e+04,3e+04)`, in bold at 1.1 times the axes font size of 24 points. A panel is 0.157 of the figure width and the panels are 0.206 apart, so at the size of `docs/images/fit.png` (12 by 6.6 inches) a panel is 135 points wide while the titles were 200 to 450 points wide. Every title ran into the next one.

## What changed

| File | Change |
| --- | --- |
| `parameterTitleText.m` (new) | Title text on two lines: `r = 1.73` and `95% CI (1.49, 2.08)`. Three significant digits; no exponent from 0.001 up to one million. |
| `fitParameterTitles.m` (new) | Gives the titles of a figure one font size that fits the panels (rule below). Also the figure's `SizeChangedFcn`. |
| `setParameterTitle.m` (new) | Writes the title in normal weight, tags it, sizes it, and sets the resize callback if the figure has none. |
| `Run_Fit_GrowthModels.m`, `Run_Forecasting_GrowthModels.m`, `plotFit_GrowthModels.m`, `plotForecast_GrowthModels.m`, `plotFit_ReproductionNumber.m`, `plotForecast_ReproductionNumber.m` | `title(cad1)` to `title(cad4)` replaced by `setParameterTitle('r',param_r(1,1:3))` and so on, placed after the line that sets the axes font size. |
| `Run_Fit_GrowthModels.m`, `Run_Forecasting_GrowthModels.m` | Semicolons after `cad1` to `cad6` and after `line1=plot(...)`, which printed to the console. The `cad` variables are still created and saved. |
| `make_readme_figures.m` | Sizes the titles again after it resizes a figure for export. |
| `tests/test_parameterTitle.m` | Ten tests (new file). |

### Font size

For each title, with `W` the width of its axes in points, `n` the number of characters of its longer line, `F` the font size of its axes and `H` the height of the figure in points:

```
min( 0.75*F,  1.25*W/(0.6*n),  0.031*H ),  not below 6 points
```

The smallest value over the four titles is used for all of them. The second term keeps a line within 1.25 axes widths if every character were 0.6 of the font size wide; digits in Helvetica and Arial are 0.556 wide and spaces, commas and parentheses are narrower. The third term keeps two lines inside the margin above the top row. The rendered width is then read from the `Extent` of each title, and the size is reduced if a title is wider than 1.15 axes widths, which covers wider fonts.

## Results (`results/titles_output.txt`)

Text was measured with Liberation Sans, which has the metrics of Helvetica and Arial. Axes positions are those of `subplot(2,4,k)` in MATLAB. Gaps are between neighbouring titles, in points; a negative gap is an overlap.

### 1. Title text

| Case | Before | After |
| --- | --- | --- |
| Quick start, `r` | `r=1.7 (95% CI:1.5,2.1)` | `r = 1.73` / `95% CI (1.49, 2.08)` |
| Quick start, `K` | `K=2.9e+04 (95% CI:2.7e+04,3e+04)` | `K = 28500` / `95% CI (27300, 29900)` |
| Rows 13–32, `r` | `r=8e+02 (95% CI:5.7e+02,1.2e+03)` | `r = 802` / `95% CI (571, 1210)` |
| Rows 13–32, `K` | `K=1.8e+04 (95% CI:1.8e+04,1.8e+04)` | `K = 18200` / `95% CI (18000, 18500)` |

With two digits the interval of `K` for rows 13–32 read as a single value. Over 20,000 random values the largest relative rounding error of the new format is 0.00497, the limit for three significant digits being 0.005.

### 2. Layout

Three sets of values: the quick start, rows 13–32, and a set with very small and very large values (`r = 0.00123`, `K = 2.85e8`).

| Figure (points) | Before: smallest gap | Before: overlapping pairs | After: smallest gap | After: font size | After: widest title, in axes widths |
| --- | --- | --- | --- | --- | --- |
| 864 × 475 (README figures) | −137 to −188 | 3 of 3 | 54 to 62 | 11.3 to 13.4 | 1.00 to 1.03 |
| 560 × 420 (default figure) | −200 to −251 | 3 of 3 | 35 to 40 | 7.3 to 8.7 | 1.00 to 1.02 |
| 1440 × 800 | −18 to −70 | 1 or 2 of 3 | 111 to 131 | 18.0 | 0.83 to 0.96 |
| 400 × 420 | −233 to −284 | 3 of 3 | 20 to 25 | 6.0 to 6.2 | 1.02 to 1.15 |

In all twelve cases the new titles are inside the figure on the right and at the top.

### 3. Wider fonts

The set with the longest titles, at the first two figure sizes.

| Font | Rule alone: widest title | After the correction by rendered width |
| --- | --- | --- |
| DejaVu Sans (about 14% wider) | 1.13 axes widths, gap 31 to 47 | not needed |
| DejaVu Sans Mono | 1.25 axes widths, gap 19 to 29 | 1.15 axes widths, gap 27 to 41, font 6.7 and 10.3 points |

### 4. Random values

Two thousand sets of four parameters with magnitudes from 1e-6 to 1e9, 5% of them negative, rendered at 864 × 475 and 560 × 420: no overlap and no title outside the figure in 4,000 layouts. Smallest gap 24.4 points, widest title 1.04 axes widths, smallest font 6.8 points.

### 5. Picture

`results/titles_before_after_mock.png` shows the row before and after at the size of the README figures. It is drawn with matplotlib from the quick-start values and is not MATLAB output.

## Limits

- **Small titles in small windows.** In the default 560 × 420 figure the titles are 7 to 9 points, next to tick labels of 24 points. Four panels of that width leave no room for more; enlarging the figure enlarges the titles, up to 18 points.
- **Other text.** Axis labels and tick labels keep the font size of `GetAdjustedFontSize`. They were not changed.
- **Numbers from matplotlib.** The gaps and font sizes above come from matplotlib, not MATLAB. MATLAB may narrow an axes to make room for its tick labels; the rule uses the width MATLAB reports. In MATLAB the layout is checked by `tests/test_parameterTitle.m`, which measures the rendered titles at three figure sizes.
- **Figures drawn before this change** and saved as `.fig` keep their old titles.

## Reproducing

```bash
pip install numpy matplotlib
python docs/validation/python/test_titles.py /tmp/titles.json
python docs/validation/python/mock_titles.py /tmp/titles_before_after_mock.png
```

The numbers above need the Liberation Sans and DejaVu fonts; with other fonts matplotlib substitutes and the widths differ.

In MATLAB, from the code folder: `results = runtests('tests/test_parameterTitle.m')`.
