"""
Regenerate the figures and every number quoted in the JASMS application note.

Usage (from the repository root):
    python paper/jasms/make_figures.py

Inputs are the result tables in paper/jasms/data/, written by cluster/sweep.slurm,
cluster/extras.slurm, the Skyline export, and cluster/procA_replicates.slurm.
Outputs: paper/jasms/figures/Figure1-3 (PNG and PDF), FigureS1, and numbers.txt.
"""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
D2, DP, OUT = HERE / "data" / "fetuin_2AB", HERE / "data" / "procainamide", HERE / "figures"
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 3, "ytick.major.size": 3,
    "lines.linewidth": 1.4, "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42})
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
SEQ = {1: "#a9cbf2", 2: "#5f9ee6", 3: "#2a78d6", 4: "#164a8a"}
MARK = ["o", "s", "^"]
ORDER = ["A2G2S1", "A2G2S2", "A3G3S1", "A3G3S2_early", "A3G3S2_late", "A3G3S3_early", "A3G3S3_late", "A3G3S4"]
LABEL = {g: g.replace("_early", " (early)").replace("_late", " (late)") for g in ORDER}
COMP = {"A2G2S1": (2, 1), "A2G2S2": (2, 2), "A3G3S1": (3, 1), "A3G3S2_early": (3, 2), "A3G3S2_late": (3, 2),
        "A3G3S3_early": (3, 3), "A3G3S3_late": (3, 3), "A3G3S4": (3, 4)}
SIALIC = {"A2G1": 0, "A2G2": 0, "A3G3": 0, "FA2G2": 0, "A2G1S1": 1, "A2G2S1": 1, "A3G3S1": 1, "FA2G2S1": 1,
          "A2G2S2": 2, "A3G3S2": 2, "FA2G2S2": 2, "A3G3S3": 3, "A3G3S4": 4}
lines = []


def say(text):
    lines.append(text)
    print(text)


def tag(ax, letter, x=-0.16):
    ax.text(x, 1.04, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", color=INK)


def save(fig, name):
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)


# ------------------------------------------------------------------ 2-AB fetuin numbers
ag = pd.read_csv(D2 / "fetuin_fld_agreement_all_settings.csv")
gl = pd.read_csv(D2 / "fetuin_glycans.csv")
pg = pd.read_csv(D2 / "fetuin_10ppm_w5_b0.1_per_glycan.csv")
me2 = pd.read_csv(D2 / "mass_errors.csv")
sk = pd.read_csv(D2 / "skyline_transition_results.csv")
add = pd.read_csv(D2 / "adduct_fld_agreement.csv")
PPM = sorted(ag.ppm.unique())

say("== 2-AB fetuin (FET-2AB-CRL)")
say(f"settings: {ag.groupby(['ppm', 'width', 'baseline_window']).ngroups}")
ok = ag[ag.ppm >= 5]
say(f">=5 ppm ({len(ok)} settings): MAD {ok.mean_abs_diff_pct_points.min():.2f}-{ok.mean_abs_diff_pct_points.max():.2f} "
    f"points, r {ok.pearson_r.min():.3f}-{ok.pearson_r.max():.3f}")
two = ag[ag.ppm == 2]
say(f"2 ppm ({len(two)} settings): MAD {two.mean_abs_diff_pct_points.min():.2f}-{two.mean_abs_diff_pct_points.max():.2f}, "
    f"species matched {two.n.min()}-{two.n.max()} of 8")
spread = gl[gl.ppm >= 5].groupby("glycan").relative_abundance_pct.agg(lambda x: x.max() - x.min())
say(f"largest per-glycan range over >=5 ppm settings: {spread.max():.2f} points ({spread.idxmax()})")
say(f"FLD delay {ag.delay_min.iloc[0] * 60:.1f} s")
d = pg.ms_pct - pg.fld_pct
say(f"default setting (10 ppm, width 5, bw 0.1): MAD {d.abs().mean():.2f}, max |d| {d.abs().max():.2f}, "
    f"LoA +/-{1.96 * d.std(ddof=1):.2f}, r {np.corrcoef(pg.ms_pct, pg.fld_pct)[0, 1]:.3f}")

s = sk.groupby("Molecule").Area.sum()
s = 100 * s / s.sum()
t = pd.DataFrame({"skyline": s, "ours": pg.set_index("glycan").ms_pct, "fld": pg.set_index("glycan").fld_pct})
table1 = []
for a, b, lab in [("ours", "skyline", "LC-MS Data Processor vs Skyline"),
                  ("ours", "fld", "LC-MS Data Processor vs fluorescence"),
                  ("skyline", "fld", "Skyline vs fluorescence")]:
    dd = t[a] - t[b]
    row = dict(comparison=lab, n=len(dd), mean_abs_diff=dd.abs().mean(), max_abs_diff=dd.abs().max(),
               loa=1.96 * dd.std(ddof=1), r=np.corrcoef(t[a], t[b])[0, 1])
    table1.append(row)
    say(f"{lab}: MAD {row['mean_abs_diff']:.2f}, max {row['max_abs_diff']:.2f}, LoA +/-{row['loa']:.2f}, r {row['r']:.3f}")
pd.DataFrame(table1).round(3).to_csv(HERE / "table1_agreement.csv", index=False)
t.round(2).to_csv(HERE / "tableS_skyline_per_glycan.csv")

m = me2[me2.adducts.isin(["2H", "3H"])].dropna(subset=["error_ppm"])
for z, g in m.groupby("charge"):
    say(f"2-AB mass error {z}+: n {len(g)}, mean {g.error_ppm.mean():+.2f}, range {g.error_ppm.min():+.2f} to {g.error_ppm.max():+.2f} ppm")

tg = pd.read_csv(D2 / "fetuin_targets.csv") if (D2 / "fetuin_targets.csv").exists() else None
for _, r in add.iterrows():
    say(f"ion set {r.ions}, {r.ppm:g} ppm: MAD {r.mean_abs_diff_pct_points:.2f}")

pg["ant"] = pg.glycan.map(lambda g: COMP[g][0])
pg["neuac"] = pg.glycan.map(lambda g: COMP[g][1])
comb = pg.groupby(["ant", "neuac"])[["ms_pct", "fld_pct"]].sum().reset_index()
comb["ratio"] = comb.ms_pct / comb.fld_pct
say("MS/FLD by antennae, NeuAc: " + "; ".join(f"{a}-ant {n} NeuAc {r:.2f}" for a, n, r in comb[["ant", "neuac", "ratio"]].values))

# ------------------------------------------------------------------ Figure 1
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw={"width_ratios": [1, 1.35]})
for k, w in enumerate(sorted(ag.width.unique())):
    st = ag[ag.width == w].groupby("ppm")["mean_abs_diff_pct_points"].agg(["min", "max", "median"]).reindex(PPM)
    x = np.array(PPM) * (1 + (k - 1) * 0.06)
    a.errorbar(x, st["median"], yerr=[st["median"] - st["min"], st["max"] - st["median"]], color=CAT[k],
               marker=MARK[k], ms=4.5, mec="white", mew=0.6, capsize=2, elinewidth=0.8, label=f"{w} points")
a.set_xscale("log")
a.set_xticks(PPM)
a.set_xticklabels([f"{p:g}" for p in PPM])
a.minorticks_off()
a.set_xlabel("Extraction window (± ppm)")
a.set_ylabel("Mean |MS − FLD| (percentage points)")
a.set_ylim(0, 8.5)
a.yaxis.grid(True, color=GRID, lw=0.5)
a.set_axisbelow(True)
a.annotate(f"{two.n.min()}–{two.n.max()} of 8 species\nmatched at ±2 ppm", xy=(2.1, 7.3), xytext=(4.5, 7.7),
           fontsize=6.5, color=INK2, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
a.legend(title="Minimum peak width", title_fontsize=7, frameon=False, loc="lower right")
tag(a, "A")

y = np.arange(len(ORDER))[::-1]
okg = gl[gl.ppm >= 5].groupby("glycan")["relative_abundance_pct"].agg(["min", "max", "median"]).reindex(ORDER)
twog = gl[gl.ppm == 2].groupby("glycan")["relative_abundance_pct"].agg(["min", "max", "median"]).reindex(ORDER)
fld = pg.set_index("glycan").reindex(ORDER)["fld_pct"]
b.hlines(y + 0.17, okg["min"], okg["max"], color=CAT[0], lw=2.2)
b.plot(okg["median"], y + 0.17, "o", color=CAT[0], ms=4.5, mec="white", mew=0.6, label="MS, ±5 to ±100 ppm (60 settings)")
b.hlines(y - 0.17, twog["min"], twog["max"], color=CAT[1], lw=2.2)
b.plot(twog["median"], y - 0.17, "s", color=CAT[1], ms=4.2, mec="white", mew=0.6, label="MS, ±2 ppm (12 settings)")
b.plot(fld, y, "D", mfc="white", mec=INK, mew=0.9, ms=4.2, label="Fluorescence")
b.set_yticks(y)
b.set_yticklabels([LABEL[g] for g in ORDER])
b.set_xlabel("Relative abundance (%)")
b.set_xlim(0, 36)
b.xaxis.grid(True, color=GRID, lw=0.5)
b.set_axisbelow(True)
b.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.18), ncol=3, handletextpad=0.3, columnspacing=1.0)
tag(b, "B", x=-0.3)
fig.tight_layout(w_pad=2.5)
save(fig, "Figure1_processing_settings")

# ------------------------------------------------------------------ Figure 2
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.1))
lim = 36
a.plot([0, lim], [0, lim], color=INK2, lw=0.8, ls="--", zorder=0)
NUM = {g: i + 1 for i, g in enumerate(ORDER)}
offsets = {"A2G2S1": (-5, 0), "A2G2S2": (6, 0), "A3G3S1": (-5, 2), "A3G3S2_early": (6, -1),
           "A3G3S2_late": (-6, 0), "A3G3S3_early": (6, 0), "A3G3S3_late": (6, -2), "A3G3S4": (6, 0)}
for g in ORDER:
    r = pg[pg.glycan == g].iloc[0]
    ant, nsia = COMP[g]
    a.plot(r.fld_pct, r.ms_pct, "o" if ant == 2 else "^", color=SEQ[nsia], ms=6.5, mec="white", mew=0.7, zorder=3)
    a.plot(r.fld_pct, t.loc[g, "skyline"], "x", color=INK2, ms=4, mew=0.8, zorder=2)
    dx, dy = offsets[g]
    a.annotate(str(NUM[g]), (r.fld_pct, r.ms_pct), xytext=(dx, dy), textcoords="offset points", fontsize=6.5,
               color=INK, va="center", ha="right" if dx < 0 else "left")
a.set_xlim(0, lim)
a.set_ylim(0, lim)
a.set_aspect("equal")
a.set_xlabel("Fluorescence relative abundance (%)")
a.set_ylabel("MS relative abundance (%)")
sw = [Line2D([], [], ls="", marker="s", color=SEQ[k], ms=5, label=f"{k}") for k in (1, 2, 3, 4)]
a.legend(handles=sw + [Line2D([], [], ls="", marker="o", color=INK2, ms=5, label="bi"),
                       Line2D([], [], ls="", marker="^", color=INK2, ms=5, label="tri"),
                       Line2D([], [], ls="", marker="x", color=INK2, ms=4, label="Skyline")],
         title="NeuAc  /  antennae", title_fontsize=6.3, ncol=3, frameon=False, loc="upper left", fontsize=6.3,
         columnspacing=0.8, handletextpad=0.2)
tag(a, "A")

b.axhline(1, color=INK2, lw=0.8, ls="--", zorder=0)
for ant, mk, col, lab in [(3, "^", CAT[0], "Triantennary (isomers combined)"), (2, "o", CAT[1], "Biantennary")]:
    sub = comb[comb.ant == ant]
    b.plot(sub.neuac, sub.ratio, color=col, marker=mk, ms=6, mec="white", mew=0.7)
    for _, r in sub.iterrows():
        b.annotate(f"{r.ratio:.2f}", (r.neuac, r.ratio), xytext=(0, 7) if ant == 2 else (0, -13),
                   textcoords="offset points", fontsize=6.3, color=INK2, ha="center")
for _, r in pg[pg.glycan.str.contains("_early|_late")].iterrows():
    early = "early" in r.glycan
    b.plot(r.neuac + (-0.22 if early else 0.22), r.ms_pct / r.fld_pct, "^", mfc="white" if early else CAT[0],
           mec=CAT[0], mew=0.8, ms=4, zorder=1)
b.set_yscale("log", base=2)
b.set_yticks([0.5, 1, 2])
b.set_yticklabels(["0.5", "1", "2"])
b.minorticks_off()
b.set_ylim(0.4, 2.8)
b.set_xticks([1, 2, 3, 4])
b.set_xlim(0.6, 4.5)
b.set_xlabel("Sialic acids (NeuAc) per glycan")
b.set_ylabel("MS share / fluorescence share")
b.yaxis.grid(True, color=GRID, lw=0.5)
b.set_axisbelow(True)
b.legend(handles=[Line2D([], [], color=CAT[0], marker="^", ms=6, mec="white", label="Triantennary (isomers combined)"),
                  Line2D([], [], color=CAT[1], marker="o", ms=6, mec="white", label="Biantennary"),
                  Line2D([], [], ls="", marker="^", mfc="white", mec=CAT[0], ms=4, label="Early-eluting isomer"),
                  Line2D([], [], ls="", marker="^", mfc=CAT[0], mec=CAT[0], ms=4, label="Late-eluting isomer")],
         frameon=False, loc="upper right", fontsize=6.3)
tag(b, "B")
fig.tight_layout(w_pad=2.5)
save(fig, "Figure2_ms_vs_fluorescence")

# ------------------------------------------------------------------ procainamide numbers
say("\n== Procainamide runs")
mp = pd.read_csv(DP / "mass_errors_all.csv").dropna(subset=["error_ppm"])
dex = pd.read_csv(DP / "dextran_ladder_DEX-ProA.csv").dropna(subset=["error_ppm"])
dex["sample"] = "DEX-ProA"
MINC = 20
for nm, g in pd.concat([mp, dex]).groupby("sample", sort=False):
    w = g[g.n_centroids >= MINC]
    say(f"{nm}: all ions n {len(g)} median {g.error_ppm.median():+.1f}; ions with >= {MINC} centroids n {len(w)}, "
        f"median {w.error_ppm.median():+.1f}, range {w.error_ppm.min():+.1f} to {w.error_ppm.max():+.1f} ppm")
wg = mp[mp.n_centroids >= MINC]
say(f"all glycan ions >= {MINC} centroids: n {len(wg)}, median {wg.error_ppm.median():+.1f}, "
    f"{(wg.error_ppm < 0).sum()} of {len(wg)} negative, {(wg.error_ppm.abs() > 10).sum()} beyond 10 ppm")
w2 = m[m.n_centroids >= 3]
say(f"2-AB ions (2+/3+): n {len(m)}, all within {m.error_ppm.min():+.2f} to {m.error_ppm.max():+.2f} ppm")

fc = pd.read_csv(DP / "fetuin_comparison.csv")
tf = pd.read_csv(DP / "transferrin_profile.csv")
for ppm in (10, 20):
    f = fc[fc.ppm == ppm].set_index("glycan")
    say(f"fetuin {ppm} ppm: " + ", ".join(f"{g} {f.loc[g, 'FET-ProA-S1']:.1f}/{f.loc[g, 'FET-ProA-S2']:.1f}" for g in f.index))
for nm in ("FET-ProA-S1", "FET-ProA-S2"):
    w = fc.pivot(index="glycan", columns="ppm", values=nm)
    dd = (w[20.0] - w[10.0]).abs()
    say(f"{nm}: largest change 10->20 ppm {dd.max():.1f} points ({dd.idxmax()})")
w = tf.pivot(index="glycan", columns="ppm", values="relative_abundance_pct")
dd = (w[20.0] - w[10.0]).abs()
say(f"TFN-ProA: largest change 10->20 ppm {dd.max():.1f} points ({dd.idxmax()})")
f20 = fc[fc.ppm == 20].set_index("glycan")
for nm in ("FET-ProA-S1", "FET-ProA-S2"):
    sia = (f20[nm] * f20.index.map(SIALIC)).sum() / 100
    asialo = f20.loc[["A2G2", "A3G3"], nm].sum()
    say(f"{nm} (20 ppm): {sia:.2f} NeuAc per glycan, asialo {asialo:.1f}%")
t20 = tf[tf.ppm == 20].set_index("glycan")
sia = (t20.relative_abundance_pct * t20.index.map(SIALIC)).sum() / 100
say(f"TFN-ProA (20 ppm): {sia:.2f} NeuAc per glycan; " + ", ".join(f"{g} {v:.1f}" for g, v in t20.relative_abundance_pct.sort_values(ascending=False).items()))

qc = pd.read_csv(DP / "qc_summary.csv").dropna(subset=["procainamide_signal"])
inc, exc = qc[qc.role == "included"], qc[qc.role == "excluded"]
say(f"QC: included signal {inc.procainamide_signal.min():.2e}-{inc.procainamide_signal.max():.2e}, ACN ratio {inc.acn_to_procainamide.min():.3f}-{inc.acn_to_procainamide.max():.3f}; "
    f"excluded signal {exc.procainamide_signal.min():.2e}-{exc.procainamide_signal.max():.2e}, ACN ratio {exc.acn_to_procainamide.min():.3f}-{exc.acn_to_procainamide.max():.3f}; "
    f"gap {inc.procainamide_signal.min() / exc.procainamide_signal.max():.1f}-fold; 2-AB signal max {qc['2AB_signal'].max()}")

# ------------------------------------------------------------------ Figure 3
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1.05, 1.15, 1]})
a, b, c = axs
runs = [("FET-2AB-CRL", m.assign(sample="FET-2AB-CRL")), ("FET-ProA-S1", mp[mp["sample"] == "FET-ProA-S1"]),
        ("FET-ProA-S2", mp[mp["sample"] == "FET-ProA-S2"]), ("TFN-ProA", mp[mp["sample"] == "TFN-ProA"]),
        ("DEX-ProA", dex)]
a.axhspan(-10, 10, color="#efeeea", zorder=0)
a.axhline(0, color=INK2, lw=0.6)
rng = np.random.default_rng(1)
for i, (nm, g) in enumerate(runs):
    good = g[g.n_centroids >= (3 if nm == "FET-2AB-CRL" else MINC)]
    col = CAT[1] if nm == "FET-2AB-CRL" else CAT[0]
    a.plot(i + rng.uniform(-0.18, 0.18, len(good)), good.error_ppm, "o", ms=3.4, color=col, mec="white", mew=0.4, zorder=3)
    a.plot([i - 0.28, i + 0.28], [good.error_ppm.median()] * 2, color=INK, lw=1.2, zorder=4)
a.set_xticks(range(len(runs)))
a.set_xticklabels([r[0] for r in runs], fontsize=6, rotation=35, ha="right", rotation_mode="anchor")
a.set_ylabel("Mass error (ppm)")
a.set_ylim(-21, 8)
a.text(4.45, 9.2, "±10 ppm window", fontsize=6, color=INK2, ha="right", va="top")
a.legend(handles=[Line2D([], [], ls="", marker="o", color=CAT[1], ms=4, label="2-AB, Agilent Q-TOF"),
                  Line2D([], [], ls="", marker="o", color=CAT[0], ms=4, label="Procainamide, SLIM IM-Q-TOF")],
         frameon=False, loc="lower left", fontsize=6, handletextpad=0.2)
tag(a, "A", x=-0.25)

GORDER = ["A2G2", "A3G3", "A2G2S1", "A3G3S1", "A2G2S2", "A3G3S2", "A3G3S3", "A3G3S4"]
yy = np.arange(len(GORDER))[::-1]
for k, (nm, off) in enumerate([("FET-ProA-S1", 0.16), ("FET-ProA-S2", -0.16)]):
    w = fc.pivot(index="glycan", columns="ppm", values=nm).reindex(GORDER)
    b.hlines(yy + off, w[10.0], w[20.0], color=GRID, lw=1.6, zorder=1)
    b.plot(w[10.0], yy + off, MARK[k], mfc="white", mec=CAT[k], mew=0.9, ms=4, zorder=2)
    b.plot(w[20.0], yy + off, MARK[k], color=CAT[k], mec="white", mew=0.5, ms=4.6, zorder=3)
b.set_yticks(yy)
b.set_yticklabels(GORDER)
b.set_xlabel("Relative abundance (%)")
b.set_xlim(0, 46)
b.xaxis.grid(True, color=GRID, lw=0.5)
b.set_axisbelow(True)
b.legend(handles=[Line2D([], [], ls="", marker="o", color=CAT[0], ms=4.5, label="FET-ProA-S1"),
                  Line2D([], [], ls="", marker="s", color=CAT[1], ms=4.5, label="FET-ProA-S2"),
                  Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, ms=4, label="±10 ppm"),
                  Line2D([], [], ls="", marker="o", color=INK2, ms=4.5, label="±20 ppm")],
         frameon=False, loc="lower right", fontsize=6, handletextpad=0.2, ncol=2, columnspacing=0.6)
tag(b, "B", x=-0.3)

for role, mk, col, lab in [("included", "o", CAT[0], "Included"), ("excluded", "s", CAT[1], "Excluded")]:
    q = qc[qc.role == role]
    c.plot(q.procainamide_signal, q.acn_to_procainamide, mk, ms=5, color=col, mec="white", mew=0.6,
           label=f"{lab} (n = {len(q)})", ls="")
c.axvline(1e7, color=INK2, lw=0.7, ls="--")
c.axhline(0.25, color=INK2, lw=0.7, ls=":")
c.set_xscale("log")
c.set_yscale("log")
c.set_xlim(3e4, 1e10)
c.set_ylim(0.02, 20)
c.set_xlabel("Procainamide-labelled signal")
c.set_ylabel("Acetonitrile-adduct / labelled signal")
c.legend(frameon=False, loc="upper right", fontsize=6, handletextpad=0.1, borderaxespad=0.1)
tag(c, "C", x=-0.3)
fig.tight_layout(w_pad=1.6)
save(fig, "Figure3_procainamide")

# ------------------------------------------------------------------ Figure S1: synthetic window benchmark
wb = pd.read_csv(D2 / "synthetic_window_benchmark.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True)
for ax, scen, letter in [(a, "mass error", "A"), (b, wb.scenario[wb.scenario != "mass error"].iloc[0], "B")]:
    sub = wb[wb.scenario == scen]
    for k, (cond, g) in enumerate(sub.groupby("condition", sort=False)):
        ax.plot(g.window_ppm, g.recovery_pct, marker=MARK[k % 3], color=CAT[k % 3], ms=4, mec="white", mew=0.5, label=cond)
    ax.axhline(100, color=INK2, lw=0.6, ls="--")
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100])
    ax.set_xticklabels(["1", "2", "5", "10", "20", "50", "100"])
    ax.minorticks_off()
    ax.set_xlabel("Extraction window (± ppm)")
    ax.legend(frameon=False, fontsize=6.5, title=scen.capitalize(), title_fontsize=6.5)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    tag(ax, letter)
a.set_ylabel("Recovered area (% of true)")
fig.tight_layout(w_pad=2)
save(fig, "FigureS1_synthetic_window_benchmark")

(HERE / "numbers.txt").write_text("\n".join(lines) + "\n")
