function make_readme_figures(rerun)
%MAKE_README_FIGURES  Regenerate the example figures shown in README.md.
%
%   make_readme_figures          runs the README examples, then exports the figures
%   make_readme_figures(false)   only redraws and exports, from the MAT files that
%                                earlier runs saved in ./output (much faster)
%
%   Run from the 'forecasting_growthmodels code' folder. The PNG files are
%   written to ../docs/images, the folder README.md links to.
%
%   Figure                     Source                                      README section
%   forecast.png               plotForecast_GrowthModels(1,1,20,4), fig 101 Quick start
%   fit.png                    plotFit_GrowthModels(1,1,20), fig 102        Quick start (fit-only)
%   forecastPerformance.png    plotForecast_GrowthModels(1,1,20,4), fig 201 Rolling windows
%   rollingForecasts.png       plotForecast_GrowthModels(1,6,20,4), fig 400 Rolling windows
%   rollingParameters.png      plotForecast_GrowthModels(1,6,20,4), the     Rolling windows
%                              unnumbered summary figure
%   Rt.png                     plotFit_ReproductionNumber(1,1,20,...), 501  Reproduction numbers
%
%   Each runner is seeded with rng(1,'twister'), as in the README. With the
%   default options the full run takes several minutes (B = 300 bootstrap
%   refits for the fit, B = 100 per window for the forecasts).
%
%   A figure that cannot be found is reported, with the figure numbers that
%   do exist, and the remaining figures are still exported. The function
%   ends with an error listing every image that was not written.

if nargin < 1 || isempty(rerun)
    rerun = true;
end

if ~exist(fullfile(pwd, 'input'), 'dir') || ~exist('Run_Forecasting_GrowthModels', 'file')
    error('make_readme_figures:Folder', ...
        'Run this function from the ''forecasting_growthmodels code'' folder.');
end

imgDir = fullfile('..', 'docs', 'images');
if ~exist(imgDir, 'dir')
    mkdir(imgDir);
end

missing = {};

% <============================================================================>
% <================== 1. Single-window forecast (Quick start) ==================>
% <============================================================================>

if rerun
    rng(1, 'twister');
    Run_Forecasting_GrowthModels(1, 1, 20, 4);
    close all
end
plotForecast_GrowthModels(1, 1, 20, 4);

missing = exportFig(findFig(101), imgDir, 'forecast.png', [12.0 6.6], missing);
missing = exportFig(findFig(201), imgDir, 'forecastPerformance.png', [8.0 6.0], missing);
close all

% <============================================================================>
% <================== 2. Fit-only analysis ======================================>
% <============================================================================>

if rerun
    rng(1, 'twister');
    Run_Fit_GrowthModels(1, 1, 20);
    close all
end
plotFit_GrowthModels(1, 1, 20);

missing = exportFig(findFig(102), imgDir, 'fit.png', [12.0 6.6], missing);
close all

% <============================================================================>
% <================== 3. Reproduction number (uses the fit above) ===============>
% <============================================================================>

plotFit_ReproductionNumber(1, 1, 20, @options_fit, @options_Rt);

missing = exportFig(findFig(501), imgDir, 'Rt.png', [9.0 7.0], missing);
close all

% <============================================================================>
% <================== 4. Rolling windows: six windows, 2-by-3 panel =============>
% <============================================================================>

if rerun
    rng(1, 'twister');
    Run_Forecasting_GrowthModels(1, 6, 20, 4);
    close all
end
plotForecast_GrowthModels(1, 6, 20, 4);

missing = exportFig(findFig(400), imgDir, 'rollingForecasts.png', [12.0 6.6], missing);
missing = exportFig(findFig('unnumbered'), imgDir, 'rollingParameters.png', [12.0 6.6], missing);
close all

% <============================================================================>

if isempty(missing)
    fprintf('All README figures written to %s\n', imgDir);
else
    error('make_readme_figures:Missing', ...
        'Not written: %s. See the messages above for the figures that were open.', ...
        strjoin(missing, ', '));
end

end


function fig = findFig(num)
%FINDFIG  The open figure with this number, or the open figure created with a
%   bare figure() call (number below 100) for num = 'unnumbered'. Returns []
%   and prints the open figure numbers when there is no unique match.

figs = findobj(groot, 'Type', 'figure');
nums = arrayfun(@(f) f.Number, figs);

if ischar(num)
    fig = figs(nums < 100);
    label = 'an unnumbered summary figure';
else
    fig = figs(nums == num);
    label = sprintf('figure %d', num);
end

if numel(fig) ~= 1
    fprintf(2, 'Expected %s; open figures: %s\n', label, mat2str(sort(nums(:))'));
    fig = [];
end

end


function missing = exportFig(fig, imgDir, name, sizeIn, missing)
%EXPORTFIG  Resize a figure (inches, width-by-height) and save it as a PNG.

if isempty(fig)
    missing{end+1} = name;
    return
end

set(fig, 'Units', 'inches');
pos = get(fig, 'Position');
set(fig, 'Position', [pos(1:2) sizeIn], 'Color', 'white');
drawnow

file = fullfile(imgDir, name);

if ~isempty(which('exportgraphics'))          % R2020a and later
    exportgraphics(fig, file, 'Resolution', 150, 'BackgroundColor', 'white');
else
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 sizeIn], 'InvertHardcopy', 'off');
    print(fig, file, '-dpng', '-r150');
end

fprintf('wrote %s\n', file);

end
