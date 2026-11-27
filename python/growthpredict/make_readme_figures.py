"""Regenerate the example figures of the README (``make_readme_figures.m``)."""

from __future__ import annotations

import os

from .errors import GrowthPredictError
from .titles import fitParameterTitles

_DEFAULT_IMAGE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs", "images")


def make_readme_figures(rerun=True, image_dir=None, *, input_dir="./input", output_dir="./output", seed=1):
    """Regenerate the example figures shown in the README.

    ``make_readme_figures()`` runs the README examples, then exports the
    figures. ``make_readme_figures(False)`` only redraws and exports, from the
    MAT files that earlier runs saved in ``output_dir`` (much faster).

    The PNG files are written to ``image_dir`` (default: ``docs/images`` inside
    the ``python`` folder of the repository).

    ======================== ====================================== ======
    Figure                   Source                                 Number
    ======================== ====================================== ======
    forecast.png             plotForecast_GrowthModels(1,1,20,4)    101
    forecastPerformance.png  plotForecast_GrowthModels(1,1,20,4)    201
    fit.png                  plotFit_GrowthModels(1,1,20)           102
    Rt.png                   plotFit_ReproductionNumber(1,1,20,...) 501
    rollingForecasts.png     plotForecast_GrowthModels(1,6,20,4)    400
    rollingParameters.png    plotForecast_GrowthModels(1,6,20,4)    300
    ======================== ====================================== ======

    Each runner is seeded with ``rng(seed)``. With the default options the full
    run takes a few minutes (B = 300 bootstrap refits for the fit, B = 100 per
    window for the forecasts). A figure that cannot be found is reported and
    the remaining figures are still exported; the function ends with an error
    listing every image that was not written. Returns the list of files written.
    """
    import matplotlib.pyplot as plt

    from .options import options_fit
    from .plotting import plotFit_GrowthModels, plotFit_ReproductionNumber, plotForecast_GrowthModels
    from .rt import options_Rt
    from .runners import Run_Fit_GrowthModels, Run_Forecasting_GrowthModels
    from .utils import rng

    imgDir = _DEFAULT_IMAGE_DIR if image_dir is None else image_dir
    os.makedirs(imgDir, exist_ok=True)

    io = dict(input_dir=input_dir, output_dir=output_dir)
    missing = []
    written = []

    def exportFig(num, name, sizeIn):
        """Resize a figure (inches, width-by-height) and save it as a PNG."""
        if num not in plt.get_fignums():
            print("Expected figure %d; open figures: %s" % (num, sorted(plt.get_fignums())))
            missing.append(name)
            return
        fig = plt.figure(num)
        fig.set_size_inches(*sizeIn)
        fig.set_facecolor("white")
        fitParameterTitles(fig)  # the titles are sized to their panels at the export size
        path = os.path.join(imgDir, name)
        fig.savefig(path, dpi=150, facecolor="white")
        written.append(path)
        print("wrote %s" % path)

    # <============================================================================>
    # <================== 1. Single-window forecast (Quick start) ==================>
    # <============================================================================>

    if rerun:
        rng(seed)
        Run_Forecasting_GrowthModels(1, 1, 20, 4, plot=False, verbose=False, **io)
    plotForecast_GrowthModels(1, 1, 20, 4, verbose=False, **io)

    exportFig(101, "forecast.png", (12.0, 6.6))
    exportFig(201, "forecastPerformance.png", (8.0, 6.0))
    plt.close("all")

    # <============================================================================>
    # <================== 2. Fit-only analysis =====================================>
    # <============================================================================>

    if rerun:
        rng(seed)
        Run_Fit_GrowthModels(1, 1, 20, plot=False, verbose=False, **io)
    plotFit_GrowthModels(1, 1, 20, verbose=False, **io)

    exportFig(102, "fit.png", (12.0, 6.6))
    plt.close("all")

    # <============================================================================>
    # <================== 3. Reproduction number (uses the fit above) ==============>
    # <============================================================================>

    plotFit_ReproductionNumber(1, 1, 20, options_fit, options_Rt, verbose=False, **io)

    exportFig(501, "Rt.png", (9.0, 7.0))
    plt.close("all")

    # <============================================================================>
    # <================== 4. Rolling windows: six windows, 2-by-3 panel ============>
    # <============================================================================>

    if rerun:
        rng(seed)
        Run_Forecasting_GrowthModels(1, 6, 20, 4, plot=False, verbose=False, **io)
    plotForecast_GrowthModels(1, 6, 20, 4, verbose=False, **io)

    exportFig(400, "rollingForecasts.png", (12.0, 6.6))
    exportFig(300, "rollingParameters.png", (12.0, 6.6))
    plt.close("all")

    # <============================================================================>

    if missing:
        raise GrowthPredictError(
            "make_readme_figures:Missing",
            "Not written: %s. See the messages above for the figures that were open." % ", ".join(missing),
        )

    print("All README figures written to %s" % imgDir)
    return written
