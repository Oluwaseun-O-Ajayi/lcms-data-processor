"""
Table-of-contents graphic (ACS: 3.25 x 1.75 in).

Panel 1: SNFG cartoon of A2G2S2 with a reducing-end label above schematic MS and fluorescence traces.
Panel 2: schematic m/z peak with +/-10 and +/-50 ppm extraction windows.
Panel 3: asialo share of FET-ProA-S1 (S1) and FET-ProA-S2 at +/-20 and +/-50 ppm (values from numbers.txt).

Usage:  python make_toc.py   (writes figures/TOC_graphic.{pdf,png,tif})
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle, RegularPolygon

plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 6, 'axes.linewidth': 0.5,
                     'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'pdf.fonttype': 42})
OUT = Path(__file__).parent / 'figures'
W, H = 3.25, 1.75
BLUE, ORANGE = '#2F6DB5', '#E0612B'
SNFG = {'GlcNAc': '#0072BC', 'Man': '#00A651', 'Gal': '#FFD400', 'NeuAc': '#A54399'}

# asialo share (%) from numbers.txt
ASIALO = {'±20 ppm': (24.6, 38.0), '±50 ppm': (35.6, 35.9)}

fig = plt.figure(figsize=(W, H))
bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, W); bg.set_ylim(0, H); bg.axis('off')


def box(x, y, w, h):
    bg.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=0.04',
                                fc='white', ec='#333333', lw=0.6))


def ax_in(x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


# title
bg.text(W / 2, 1.62, 'Sweep the extraction window before reporting a glycan profile',
        ha='center', va='center', fontsize=7, weight='bold')

px = [0.06, 1.15, 2.24]
pw, py, ph = 0.95, 0.08, 1.38
for x in px:
    box(x, py, pw, ph)
for x0, x1 in [(px[0] + pw, px[1]), (px[1] + pw, px[2])]:
    bg.add_patch(FancyArrowPatch((x0 + 0.02, py + ph / 2), (x1 - 0.02, py + ph / 2),
                                 arrowstyle='-|>', mutation_scale=8, lw=0.6, color='#333333'))
heads = ['MS vs fluorescence', 'Window sweep', 'Profile vs window']
for x, t in zip(px, heads):
    bg.text(x + pw / 2, py + ph - 0.08, t, ha='center', va='center', fontsize=6, weight='bold')

# ---- panel 1: glycan cartoon
s = 0.055   # symbol size (in)
gy = 1.20
lab_x = px[0] + 0.11
bg.add_patch(FancyBboxPatch((lab_x - 0.07, gy - 0.032), 0.14, 0.064, boxstyle='round,pad=0,rounding_size=0.015',
                            fc='#BBBBBB', ec='#555555', lw=0.4))
bg.text(lab_x, gy, 'Label', ha='center', va='center', fontsize=3.7)
core = [lab_x + 0.14, lab_x + 0.23]
man_c = lab_x + 0.32
arms = [gy + 0.07, gy - 0.07]
lines = [(lab_x + 0.07, gy, core[0], gy), (core[0], gy, core[1], gy), (core[1], gy, man_c, gy)]
for ay in arms:
    xs = [man_c + 0.09, man_c + 0.18, man_c + 0.27, man_c + 0.36]
    lines += [(man_c, gy, xs[0], ay)] + [(xs[i], ay, xs[i + 1], ay) for i in range(3)]
for x0, y0, x1, y1 in lines:
    bg.plot([x0, x1], [y0, y1], color='#333333', lw=0.5, zorder=1)


def sq(x, y):
    bg.add_patch(Rectangle((x - s / 2, y - s / 2), s, s, fc=SNFG['GlcNAc'], ec='black', lw=0.35, zorder=2))


def ci(x, y, c):
    bg.add_patch(Circle((x, y), s / 2, fc=SNFG[c], ec='black', lw=0.35, zorder=2))


def di(x, y):
    bg.add_patch(RegularPolygon((x, y), 4, radius=s * 0.62, orientation=0, fc=SNFG['NeuAc'], ec='black',
                                lw=0.35, zorder=2))


sq(core[0], gy); sq(core[1], gy); ci(man_c, gy, 'Man')
for ay in arms:
    ci(man_c + 0.09, ay, 'Man'); sq(man_c + 0.18, ay); ci(man_c + 0.27, ay, 'Gal'); di(man_c + 0.36, ay)

# traces
a1 = ax_in(px[0] + 0.14, py + 0.14, pw - 0.20, 0.78)
t = np.linspace(0, 10, 600)
ms = 1.0 * np.exp(-((t - 5) / 0.35) ** 2) + 0.12 + 0.01 * np.sin(3 * t)
fl = 0.72 * np.exp(-((t - 5.05) / 0.38) ** 2) + 0.06
a1.plot(t, ms, color=BLUE, lw=0.8, label='MS (XIC)')
a1.plot(t, fl, color=ORANGE, lw=0.8, label='Fluorescence')
a1.set_xticks([]); a1.set_yticks([]); a1.set_ylim(0, 1.45)
a1.set_xlabel('Time', labelpad=1); a1.set_ylabel('Intensity', labelpad=1)
a1.legend(loc='upper right', fontsize=4.5, frameon=False, handlelength=1.0, borderaxespad=0.1)
for sp in ('top', 'right'):
    a1.spines[sp].set_visible(False)

# ---- panel 2: windows on an m/z peak
a2 = ax_in(px[1] + 0.14, py + 0.14, pw - 0.20, 1.02)
mz = np.linspace(-80, 80, 600)
a2.axvspan(-50, -10, color=ORANGE, alpha=0.25, lw=0)
a2.axvspan(10, 50, color=ORANGE, alpha=0.25, lw=0)
a2.axvspan(-10, 10, color=BLUE, alpha=0.35, lw=0)
a2.plot(mz, np.exp(-((mz + 6) / 22) ** 2), color='black', lw=0.8)
a2.text(0, 1.12, '±10 ppm', ha='center', fontsize=5, color=BLUE, weight='bold')
a2.text(48, 0.80, '±50\nppm', ha='right', va='center', fontsize=5, color=ORANGE, weight='bold', linespacing=0.95)
a2.set_ylim(0, 1.25); a2.set_xlim(-80, 80)
a2.set_xticks([]); a2.set_yticks([])
a2.set_xlabel('m/z', labelpad=1, style='italic'); a2.set_ylabel('Intensity', labelpad=1)
for sp in ('top', 'right'):
    a2.spines[sp].set_visible(False)

# ---- panel 3: asialo share, real values
a3 = ax_in(px[2] + 0.22, py + 0.17, pw - 0.28, 0.96)
xpos = []
for i, (k, (s1, s2)) in enumerate(ASIALO.items()):
    x0 = i * 1.0
    a3.bar([x0 - 0.18, x0 + 0.18], [s1, s2], width=0.34, color=[BLUE, ORANGE], ec='none')
    xpos.append(x0)
a3.set_xticks(xpos); a3.set_xticklabels(list(ASIALO), fontsize=5)
a3.tick_params(axis='x', length=0, pad=1.5)
a3.set_ylim(0, 45); a3.set_yticks([0, 20, 40]); a3.tick_params(axis='y', labelsize=5, length=1.5, pad=1)
a3.set_ylabel('Asialo share (%)', labelpad=1, fontsize=5.5)
for sp in ('top', 'right'):
    a3.spines[sp].set_visible(False)
a3.bar([0], [0], color=BLUE, label='S1'); a3.bar([0], [0], color=ORANGE, label='S2')
a3.legend(loc='upper left', fontsize=4.5, frameon=False, handlelength=0.8, borderaxespad=0.1, ncol=2,
          columnspacing=0.6, bbox_to_anchor=(0.0, 1.10))

OUT.mkdir(exist_ok=True)
fig.savefig(OUT / 'TOC_graphic.pdf')
fig.savefig(OUT / 'TOC_graphic.png', dpi=600)
fig.savefig(OUT / 'TOC_graphic.tif', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
print('wrote TOC_graphic.{pdf,png,tif}')
