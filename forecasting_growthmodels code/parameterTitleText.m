function lines = parameterTitleText(name, values)
%PARAMETERTITLETEXT  Two-line title text for a parameter histogram.
%
%   lines = parameterTitleText(name, values)
%
%   name    parameter name, e.g. 'r'
%   values  [estimate lower upper], the estimate and the limits of its 95% CI
%
%   lines   1-by-2 cell, e.g. {'r = 1.73', '95% CI (1.49, 2.08)'}
%
%   Numbers are written with three significant digits. Values with a
%   magnitude from 0.001 up to (not including) 1e6 are written without an
%   exponent and without trailing zeros (28547 -> '28500', 0.8526 ->
%   '0.853', 1 -> '1'); others as '2.85e8' or '3.20e-5'.
%
%   The estimate goes on the first line and the interval on the second so
%   that each line is short enough to fit above one panel of a row of four.
%   See also setParameterTitle, fitParameterTitles.

if ~(ischar(name) || (isstring(name) && isscalar(name)))
    error('GrowthPredict:parameterTitleText:name', ...
        'name must be a character vector, e.g. ''r''.');
end

if ~isnumeric(values) || numel(values) ~= 3
    error('GrowthPredict:parameterTitleText:values', ...
        'values must be [estimate lower upper].');
end

values = double(values(:)');

lines = {sprintf('%s = %s', char(name), formatNumber(values(1))), ...
    sprintf('95%% CI (%s, %s)', formatNumber(values(2)), formatNumber(values(3)))};

end


function s = formatNumber(x)
%FORMATNUMBER  Three significant digits; no exponent between 1e-3 and 1e6.

if isnan(x)
    s = 'NaN';
    return
end

if isinf(x)
    if x > 0
        s = 'Inf';
    else
        s = '-Inf';
    end
    return
end

if x == 0
    s = '0';
    return
end

y = str2double(sprintf('%.3g', x));      % x rounded to three significant digits
ay = abs(y);

if ay >= 1e-3 && ay < 1e6

    decimals = max(0, 2 - floor(log10(ay)));
    s = sprintf('%.*f', decimals, y);

    if any(s == '.')
        s = regexprep(s, '0+$', '');     % trailing zeros
        s = regexprep(s, '\.$', '');     % and a decimal point left at the end
    end

else

    s = sprintf('%.2e', y);              % e.g. 2.85e+08
    s = regexprep(s, 'e([+-])0*(\d)', 'e$1$2');
    s = strrep(s, 'e+', 'e');            % 2.85e8

end

end
