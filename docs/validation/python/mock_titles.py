"""Mock of the top row of fit.png (12 x 6.6 in, 24-pt axes font) drawn with matplotlib, before and after.
Not MATLAB output: it shows the title layout only."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from title_port import title_lines, old_title, rule_font_size
LEFT = [0.1300, 0.3361, 0.5422, 0.7484]; WIDTH = 0.1566; W, H = 864, 475.2
S = dict(r=(1.7302, 1.49, 2.08), p=(0.8526, 0.828, 0.874), a=(1, 1, 1), K=(28547, 27295, 29915))
rng = np.random.default_rng(0)
fig = plt.figure(figsize=(W/72, 2*H/72*0.62), dpi=72)
for row, new in enumerate((False, True)):
    lines = [title_lines(k, S[k]) for k in 'rpaK']; fs = rule_font_size([WIDTH*W]*4, [24]*4, lines, H)
    for i, k in enumerate('rpaK'):
        ax = fig.add_axes([LEFT[i], 0.60 - 0.50*row, WIDTH, 0.25])
        est, lo, hi = S[k]; x = rng.normal(est, max((hi-lo)/4, 1e-9), 300) if hi > lo else np.full(300, est)
        ax.hist(x, color='#3b6fb6', edgecolor='white'); ax.axvline(lo, color='r', ls='--', lw=2); ax.axvline(hi, color='r', ls='--', lw=2); ax.axvline(est, color='k', ls='--', lw=2)
        ax.set_xlabel(k, fontsize=20, family='Liberation Sans'); ax.tick_params(labelsize=13)
        if i == 0: ax.set_ylabel('Frequency', fontsize=20, family='Liberation Sans')
        if new: ax.set_title("\n".join(lines[i]), fontproperties=FontProperties(family='Liberation Sans', size=fs))
        else:   ax.set_title(old_title(k, S[k]), fontproperties=FontProperties(family='Liberation Sans', weight='bold', size=26.4))
    fig.text(0.02, 0.97 - 0.50*row, 'Before' if not new else 'After', fontsize=15, family='Liberation Sans', color='#555555', va='top')
fig.savefig(sys.argv[1], dpi=100, facecolor='white')
