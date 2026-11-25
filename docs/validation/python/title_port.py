"""Python port of parameterTitleText.m and of the font-size rule in fitParameterTitles.m."""
import math, re

CAP, ROOM, CHAR_WIDTH, HEIGHT_FRACTION, MAX_FRACTION, MIN_FONT = 0.75, 1.25, 0.6, 0.031, 1.15, 6.0


def format_number(x):
    """Three significant digits, no exponent between 1e-3 and 1e6, trailing zeros removed."""
    if isinstance(x, float) and math.isnan(x):
        return 'NaN'
    if isinstance(x, float) and math.isinf(x):
        return 'Inf' if x > 0 else '-Inf'
    if x == 0:
        return '0'
    y = float(f"{x:.3g}")
    ay = abs(y)
    if 1e-3 <= ay < 1e6:
        decimals = max(0, 2 - math.floor(math.log10(ay)))
        s = f"{y:.{decimals}f}"
        if '.' in s:
            s = s.rstrip('0').rstrip('.')
        return s
    s = f"{y:.2e}"
    return re.sub(r'e([+-])0*(\d)', r'e\1\2', s).replace('e+', 'e')


def title_lines(name, values):
    est, lo, hi = values
    return [f"{name} = {format_number(est)}", f"95% CI ({format_number(lo)}, {format_number(hi)})"]


def old_title(name, values):
    """Runners before the change: strcat(name,'=',num2str(est,2),' (95% CI:',num2str(lo,2),',',num2str(hi,2),')')."""
    n2 = lambda v: re.sub(r'e([+-])(\d)$', r'e\g<1>0\2', f"{v:.2g}")       # MATLAB prints a two-digit exponent
    return f"{name}={n2(values[0])} (95% CI:{n2(values[1])},{n2(values[2])})"


def rule_font_size(axes_widths_pt, axes_fonts_pt, all_lines, fig_height_pt):
    """One font size for all titles of a figure: the smallest of the per-title values, not below MIN_FONT."""
    sizes = []
    for w, f, lines in zip(axes_widths_pt, axes_fonts_pt, all_lines):
        nchar = max(1, max(len(s) for s in lines))
        sizes.append(min(CAP * f, ROOM * w / (CHAR_WIDTH * nchar), HEIGHT_FRACTION * fig_height_pt))
    return max(MIN_FONT, min(sizes))


def corrected_font_size(font_size, widest_fraction):
    """Second step: reduce the size if a rendered title is wider than MAX_FRACTION axes widths."""
    if math.isfinite(widest_fraction) and widest_fraction > MAX_FRACTION:
        return max(MIN_FONT, font_size * MAX_FRACTION / widest_fraction)
    return font_size
