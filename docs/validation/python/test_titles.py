"""Panel titles of the parameter histograms: one-line bold title (before) vs two-line fitted title (after).

Text is measured with Liberation Sans, which has the metrics of Helvetica and Arial (MATLAB's default
fonts on macOS and Windows), and with DejaVu Sans, which is about 14% wider, to check the correction
by rendered width. Axes positions are those of MATLAB's subplot(2,4,k).

usage: python3 test_titles.py results.json [folder for mock figures]
"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from title_port import format_number, title_lines, old_title, rule_font_size, corrected_font_size, MIN_FONT
OUT = {}
LEFT = [0.1300, 0.3361, 0.5422, 0.7484]; WIDTH = 0.1566; BOTTOM = 0.5838; HEIGHT = 0.3412     # MATLAB subplot(2,4,k), k = 1..4
AXFONT = 24.0                                                                                # GetAdjustedFontSize
fails = 0
def check(ok, msg):
    global fails; fails += (not ok); print(f"  [{'PASS' if ok else 'FAIL'}] {msg}")

print("=== 1. Number format (three significant digits)")
cases = {1.7302: '1.73', 0.8526: '0.853', 1: '1', 28547.1: '28500', 801.8: '802', 0.17743: '0.177', 18153.2: '18200', 0.0012345: '0.00123', 999.6: '1000', 9.996: '10',
         123456.7: '123000', 1234567.0: '1.23e6', 2.85e8: '2.85e8', 3.2e-5: '3.20e-5', 0: '0', -0.52: '-0.52', 10: '10', 0.5: '0.5', 100: '100', 999999.5: '1.00e6', 0.001: '0.001', 0.00099949: '9.99e-4'}
for c, e in cases.items(): check(format_number(c) == e, f"{c!r} -> {format_number(c)}")
check(format_number(float('nan')) == 'NaN' and format_number(float('inf')) == 'Inf' and format_number(float('-inf')) == '-Inf', "NaN, Inf, -Inf")
rng = np.random.default_rng(1); worst = 0.0
for _ in range(20000):
    x = 10**rng.uniform(-8, 10) * rng.choice([-1, 1]); worst = max(worst, abs(float(format_number(x)) - x)/abs(x))
check(worst <= 0.005001, f"20,000 random values: largest relative rounding error {worst:.5f} (three significant digits allow 0.005)")
OUT['numbers'] = {repr(c): format_number(c) for c in cases}

SETS = {'quick start (rows 1-20)': dict(r=(1.7302, 1.49, 2.08), p=(0.8526, 0.828, 0.874), a=(1, 1, 1), K=(28547, 27295, 29915)),
        'rows 13-32':              dict(r=(801.8, 571, 1210), p=(0.178, 0.128, 0.223), a=(1, 1, 1), K=(18153, 18026, 18470)),
        'large and small values':  dict(r=(0.00123, 0.00101, 0.00187), p=(0.999, 0.951, 1), a=(4.07, 0.0123, 9.87), K=(2.85e8, 1.02e8, 9.91e8))}
print("\n=== 2. Title text")
for nm, S in SETS.items():
    for k, v in S.items():
        print(f"  {nm:<26} before: {old_title(k, v):<34} after: {' | '.join(title_lines(k, v))}")

def measure(fig_w_pt, fig_h_pt, S, new, family='Liberation Sans', correct=True):
    """Render the four titles; return the figure, the edges (pt) of each title and the font size."""
    fig = plt.figure(figsize=(fig_w_pt/72, fig_h_pt/72), dpi=72); r = fig.canvas.get_renderer()
    axs = []
    for i in range(4):
        ax = fig.add_axes([LEFT[i], BOTTOM, WIDTH, HEIGHT]); ax.set_xticks([]); ax.set_yticks([]); axs.append(ax)
    def draw(fs, weight, texts):
        ts = [ax.set_title(t, fontproperties=FontProperties(family=family, weight=weight, size=fs)) for ax, t in zip(axs, texts)]
        fig.canvas.draw(); return [t.get_window_extent(r) for t in ts]
    if new:
        lines = [title_lines(k, S[k]) for k in 'rpaK']
        fs = rule_font_size([WIDTH*fig_w_pt]*4, [AXFONT]*4, lines, fig_h_pt)
        bbs = draw(fs, 'normal', ["\n".join(l) for l in lines])
        if correct:
            fs2 = corrected_font_size(fs, max(b.width for b in bbs)/(WIDTH*fig_w_pt))
            if fs2 != fs: fs = fs2; bbs = draw(fs, 'normal', ["\n".join(l) for l in lines])
    else:
        fs = AXFONT*1.1; bbs = draw(fs, 'bold', [old_title(k, S[k]) for k in 'rpaK'])
    return fig, bbs, fs

def report(bbs, fig_w, fig_h):
    gaps = [bbs[i+1].x0 - bbs[i].x1 for i in range(3)]            # negative = overlap with the next title
    return min(gaps), fig_w - bbs[3].x1, sum(1 for g in gaps if g < 0), fig_h - max(b.y1 for b in bbs)

print("\n=== 3. Layout with Helvetica/Arial metrics: smallest gap between neighbouring titles (points; negative = overlap)")
FIGS = {'README export 864 x 475': (864, 475.2), 'default figure 560 x 420': (560, 420), 'large 1440 x 800': (1440, 800), 'narrow 400 x 420': (400, 420)}
for fname, (W, H) in FIGS.items():
    for nm, S in SETS.items():
        fo, bo, so = measure(W, H, S, False); go, ro, no, to = report(bo, W, H)
        fn, bn, sn = measure(W, H, S, True); gn, rn, nn, tn = report(bn, W, H)
        widest = max(b.width for b in bn)/(WIDTH*W)
        print(f"  {fname:<25} {nm:<24} before: gap {go:7.1f}, overlapping pairs {no}, {so:.1f} pt bold | after: gap {gn:6.1f}, right margin {rn:6.1f}, top margin {tn:5.1f}, font {sn:4.1f} pt, widest {widest:.2f} axes widths")
        check(gn > 0 and rn > 0 and tn > 0 and nn == 0, "no overlap, inside the figure")
        OUT[f"{fname}|{nm}"] = dict(before_gap=go, before_pairs=no, after_gap=gn, after_right=rn, after_top=tn, after_pairs=nn, font=sn, widest=widest)
        if fname.startswith('README') and nm.startswith('quick') and len(sys.argv) > 2:
            fo.savefig(os.path.join(sys.argv[2], 'titles_before_mock.png'), dpi=110, facecolor='white')
            fn.savefig(os.path.join(sys.argv[2], 'titles_after_mock.png'), dpi=110, facecolor='white')
        plt.close('all')
check(all(v['before_pairs'] >= 1 for k, v in OUT.items() if '|' in k), "the one-line bold titles overlap in every case above")

print("\n=== 4. Wider fonts: rule alone, then with the correction by rendered width")
for family in ('DejaVu Sans', 'DejaVu Sans Mono'):
    for fname, (W, H) in list(FIGS.items())[:2]:
        S = SETS['large and small values']
        _, b0, s0 = measure(W, H, S, True, family, correct=False); g0 = report(b0, W, H)[0]; w0 = max(b.width for b in b0)/(WIDTH*W)
        _, b1, s1 = measure(W, H, S, True, family, correct=True);  g1 = report(b1, W, H)[0]; w1 = max(b.width for b in b1)/(WIDTH*W)
        print(f"  {family:<17} {fname:<25} rule: font {s0:.1f} pt, widest {w0:.2f} axes widths, gap {g0:.1f} | corrected: font {s1:.1f} pt, widest {w1:.2f}, gap {g1:.1f}")
        check(g1 > 0 and w1 <= 1.15 + 0.01 and s1 <= s0, "titles are at most 1.15 axes widths and do not overlap")
        OUT[f"{family}|{fname}"] = dict(rule_font=s0, rule_widest=w0, rule_gap=g0, corrected_font=s1, corrected_widest=w1, corrected_gap=g1)
        plt.close('all')
check(OUT['DejaVu Sans Mono|README export 864 x 475']['corrected_font'] < OUT['DejaVu Sans Mono|README export 864 x 475']['rule_font'], "the correction is used for the monospaced font")

print("\n=== 5. Random parameter values (2,000 sets, magnitudes 1e-6 to 1e9, 5% negative), rendered at two figure sizes")
rng = np.random.default_rng(3); bad = 0; worst = 1e9; widest = 0; fmin = 1e9
fig_cache = {}
for it in range(2000):
    S = {}
    for k in 'rpaK':
        est = 10**rng.uniform(-6, 9) * (1 if rng.random() > 0.05 else -1); a_, b_ = est*rng.uniform(0.3, 1), est*rng.uniform(1, 3)
        S[k] = (est, min(a_, b_), max(a_, b_))
    for (W, H) in ((864, 475.2), (560, 420)):
        _, bn, sn = measure(W, H, S, True); g, r_, n, t = report(bn, W, H)
        bad += (g <= 0) or (r_ <= 0) or (t <= 0); worst = min(worst, g); widest = max(widest, max(b.width for b in bn)/(WIDTH*W)); fmin = min(fmin, sn)
        plt.close('all')
print(f"  overlapping or clipped layouts: {bad} of 4000 | smallest gap {worst:.1f} pt | widest title {widest:.2f} axes widths (panels are 1.32 axes widths apart) | smallest font {fmin:.1f} pt")
check(bad == 0, "no overlap in any random layout")
OUT['random'] = dict(bad=int(bad), smallest_gap=float(worst), widest=float(widest), smallest_font=float(fmin))
json.dump(OUT, open(sys.argv[1], 'w'), indent=1)
print(f"\n{'ALL CHECKS PASSED' if fails == 0 else str(fails) + ' CHECKS FAILED'}")
sys.exit(1 if fails else 0)
