"""
Regenerate the figures and every number quoted in the released N-glycan validation paper.

Usage (from the repository root):
    python paper/glycan-validation/make_figures.py

Inputs are the result tables in paper/glycan-validation/data/, written by cluster/sweep.slurm,
cluster/extras.slurm, the Skyline export, cluster/procA_replicates.slurm (fetuin, procainamide),
combine_replicates.py, transfer_test.py and lcms-compare fld (RapiFluor-MS data set PXD078897).
Outputs: paper/glycan-validation/figures/Figure1-3 and FigureS1-S4 (PNG and PDF), TOC_graphic, and numbers.txt.
Figure S4 (isomer pairs) needs the traces written by xic_pair.py (cluster/xic_pairs.slurm).
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
DT, DM = HERE / "data" / "trastuzumab_RFMS", HERE / "data" / "ms119_transfer"
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
fig1, (a, b, c1) = plt.subplots(1, 3, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1, 1.25, 1]})
fig = fig1
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
a.set_ylabel("Mean |MS − FLD| (points)")
a.set_ylim(0, 8.5)
a.yaxis.grid(True, color=GRID, lw=0.5)
a.set_axisbelow(True)
a.annotate(f"{two.n.min()}–{two.n.max()} of 8 species\nmatched at ±2 ppm", xy=(2.1, 7.3), xytext=(4.5, 7.7),
           fontsize=6.5, color=INK2, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
a.legend(title="Minimum peak width", title_fontsize=7, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.86))
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
b.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.35, -0.2), ncol=2, handletextpad=0.3, columnspacing=1.0)
tag(b, "B", x=-0.42)
a1, b1 = a, b

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
save(fig, "FigureS3_fetuin_ms_vs_fluorescence")

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
WIN = sorted(fc.ppm.unique())
for nm in ("FET-ProA-S1", "FET-ProA-S2"):
    w = fc.pivot(index="glycan", columns="ppm", values=nm)
    steps = ", ".join(f"{a:g}->{b:g} {(w[b] - w[a]).abs().max():.1f} ({(w[b] - w[a]).abs().idxmax()})" for a, b in zip(WIN, WIN[1:]))
    say(f"{nm}: largest change between consecutive windows: {steps}; 10->50 {(w[WIN[-1]] - w[WIN[0]]).abs().max():.1f}")
    for p_ in WIN:
        sia = (w[p_] * w.index.map(SIALIC)).sum() / 100
        say(f"{nm} ±{p_:g} ppm: {sia:.2f} NeuAc per glycan, asialo {w.loc[['A2G2', 'A3G3'], p_].sum():.1f}%, "
            + ", ".join(f"{g} {w.loc[g, p_]:.1f}" for g in ["A2G2", "A3G3", "A3G3S1", "A3G3S2"]))
dd = fc.pivot(index="glycan", columns="ppm", values="difference_pct_points").abs()
say("largest |S2 - S1| by window: " + ", ".join(f"±{p_:g} {dd[p_].max():.1f} ({dd[p_].idxmax()})" for p_ in WIN))
t = tf.pivot(index="glycan", columns="ppm", values="relative_abundance_pct")
steps = ", ".join(f"{a:g}->{b:g} {(t[b] - t[a]).abs().max():.1f} ({(t[b] - t[a]).abs().idxmax()})" for a, b in zip(WIN, WIN[1:]))
say(f"TFN-ProA: largest change between consecutive windows: {steps}")
for g in ["A2G2S1", "A2G2", "A2G2S2"]:
    say(f"TFN-ProA {g}: " + ", ".join(f"±{p_:g} {t.loc[g, p_]:.1f}" for p_ in WIN))
tg = pd.read_csv(DP / "TFN-ProA_targets.csv")
a2 = tg[tg.target == "A2G2_2+"].set_index("ppm").peak_area
say("TFN-ProA A2G2 2+ area by window: " + ", ".join(f"±{p_:g} {a2[p_]:.3g}" for p_ in WIN))

qc = pd.read_csv(DP / "qc_summary.csv").dropna(subset=["procainamide_signal"])
inc, exc = qc[qc.role == "included"], qc[qc.role == "excluded"]
say(f"QC: included signal {inc.procainamide_signal.min():.2e}-{inc.procainamide_signal.max():.2e}, ACN ratio {inc.acn_to_procainamide.min():.3f}-{inc.acn_to_procainamide.max():.3f}; "
    f"excluded signal {exc.procainamide_signal.min():.2e}-{exc.procainamide_signal.max():.2e}, ACN ratio {exc.acn_to_procainamide.min():.3f}-{exc.acn_to_procainamide.max():.3f}; "
    f"gap {inc.procainamide_signal.min() / exc.procainamide_signal.max():.1f}-fold; 2-AB signal max {qc['2AB_signal'].max()}")

# ------------------------------------------------------------------ Figure 3
fig = plt.figure(figsize=(7.0, 5.0))
gs = fig.add_gridspec(2, 5, height_ratios=[1.15, 1], hspace=0.5, wspace=0.75)
a = fig.add_subplot(gs[0, 0:3])
c = fig.add_subplot(gs[0, 3:5])
c.set_position(c.get_position())
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
a.set_xticklabels([r[0] for r in runs], fontsize=6.3)
a.set_ylabel("Mass error (ppm)")
a.set_ylim(-21, 8)
a.text(4.45, 9.2, "±10 ppm window", fontsize=6, color=INK2, ha="right", va="top")
a.legend(handles=[Line2D([], [], ls="", marker="o", color=CAT[1], ms=4, label="2-AB, Agilent Q-TOF"),
                  Line2D([], [], ls="", marker="o", color=CAT[0], ms=4, label="Procainamide, SLIM IM-Q-TOF")],
         frameon=False, loc="lower left", fontsize=6, handletextpad=0.2)
tag(a, "A", x=-0.12)

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
c.set_xlabel("Procainamide-labeled signal")
c.set_ylabel("Acetonitrile-adduct ratio")
c.legend(frameon=False, loc="upper right", fontsize=6, handletextpad=0.1, borderaxespad=0.1)
tag(c, "B", x=-0.28)

panels = [("A2G2", "fet"), ("A3G3", "fet"), ("A3G3S1", "fet"), ("A3G3S2", "fet"), ("A2G2S1", "tfn")]
for k, (g, kind) in enumerate(panels):
    ax = fig.add_subplot(gs[1, k])
    if kind == "fet":
        for nm, col, mk, ls in [("FET-ProA-S1", CAT[0], "o", "-"), ("FET-ProA-S2", CAT[1], "s", "--")]:
            w = fc.pivot(index="glycan", columns="ppm", values=nm).loc[g]
            ax.plot(WIN, w[WIN].values, ls=ls, marker=mk, color=col, ms=3.5, mec="white", mew=0.4, lw=1.1)
        ax.set_title(f"Fetuin {g}", fontsize=7)
    else:
        ax.plot(WIN, t.loc["A2G2S1", WIN].values, "-^", color=CAT[2], ms=3.5, mec="white", mew=0.4, lw=1.1)
        ax.set_title("Transferrin A2G2S1", fontsize=7)
    ax.set_xticks(WIN)
    ax.set_xticklabels([f"{p_:g}" for p_ in WIN], fontsize=6)
    ax.tick_params(axis="y", labelsize=6)
    lo, hi = ax.get_ylim()
    ax.set_ylim(0, max(hi * 1.1, 5))
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    if k == 0:
        ax.set_ylabel("Relative abundance (%)")
        tag(ax, "C", x=-0.6)
    if k == 2:
        ax.set_xlabel("Extraction window (± ppm)")
fig.legend(handles=[Line2D([], [], ls="-", marker="o", color=CAT[0], ms=3.5, label="FET-ProA-S1"),
                    Line2D([], [], ls="--", marker="s", color=CAT[1], ms=3.5, label="FET-ProA-S2"),
                    Line2D([], [], ls="-", marker="^", color=CAT[2], ms=3.5, label="TFN-ProA")],
           frameon=False, loc="lower center", ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, -0.01))
fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.14)
save(fig, "FigureS2_procainamide")

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


# ------------------------------------------------------------------ RapiFluor-MS data set (PXD078897)
GLY = ["Man5", "G0", "G0F", "G1F_a", "G1F_b", "G2F"]  # in order of mass
GLAB = {"Man5": "Man5", "G0": "G0", "G0F": "G0F", "G1F_a": "G1F (early)", "G1F_b": "G1F (late)", "G2F": "G2F", "G1": "G1a"}
DAYM = {1: "o", 2: "s", 3: "^"}
pr = pd.read_csv(DT / "per_run_per_glycan.csv")
say("\n== RapiFluor-MS, MS154 trastuzumab (9 runs, 3 days)")
dd = pr.ms_pct - pr.fld_pct
say(f"pooled: n {len(dd)}, MAD {dd.abs().mean():.3f}, LoA +/-{1.96 * dd.std(ddof=1):.3f}, "
    f"r {np.corrcoef(pr.ms_pct, pr.fld_pct)[0, 1]:.4f}")
pr["rf"] = pr.ms_pct / pr.fld_pct
rf_all = pr.groupby("glycan").rf.mean().reindex(GLY)
rf_day = pr.groupby(["glycan", "day"]).rf.mean().unstack().reindex(GLY)
say("response factor, mean of 9 runs: " + ", ".join(f"{g} {v:.3f}" for g, v in rf_all.items()))
say("response factor range over days: " + ", ".join(f"{g} {r.min():.3f}-{r.max():.3f}" for g, r in rf_day.iterrows()))
bias = []
for g, x in pr.groupby("glycan"):
    half = __import__("scipy.stats", fromlist=["t"]).t.ppf(0.975, len(x) - 1) * (x.ms_pct - x.fld_pct).std(ddof=1) / np.sqrt(len(x))
    m_ = (x.ms_pct - x.fld_pct).mean()
    bias.append({"glycan": g, "bias": m_, "lo": m_ - half, "hi": m_ + half})
bias = pd.DataFrame(bias).set_index("glycan").reindex(GLY)
say("bias [95% CI]: " + ", ".join(f"{g} {r.bias:+.2f} [{r.lo:.2f}, {r.hi:.2f}]" for g, r in bias.iterrows()))
sens = pd.read_csv(DT / "settings_sensitivity.csv")
bys = sens.groupby(["ppm", "width", "baseline_window"]).mean_abs_diff.mean()
say(f"settings: mean MAD over runs {bys.min():.3f}-{bys.max():.3f} across {len(bys)} settings")
cal = pd.read_csv(DT / "calibration_cross_day.csv")
cd = cal.groupby("train_day")[["mad_before", "mad_after", "loa_halfwidth_before", "loa_halfwidth_after"]].mean()
say(f"cross-day calibration: MAD {cal.mad_before.mean():.3f} -> {cal.mad_after.mean():.3f}; by training day after "
    f"{cd.mad_after.min():.3f}-{cd.mad_after.max():.3f}; LoA half-width {cd.loa_halfwidth_before.min():.2f}-"
    f"{cd.loa_halfwidth_before.max():.2f} -> {cd.loa_halfwidth_after.min():.2f}-{cd.loa_halfwidth_after.max():.2f}")
full = pd.read_csv(DM / "full_profile.csv")
trf = pd.read_csv(DM / "transfer.csv")
wpr = pd.read_csv(DM / "within_protein.csv")
PROT = {"Trast": "Trastuzumab", "Afli": "Aflibercept", "Cetux": "Cetuximab"}
for p_, g in full.groupby("protein"):
    say(f"MS119 {p_}: n {int(g.n.iloc[0])}, MAD {g.mean_abs_diff.mean():.3f}, LoA +/-{g.loa_halfwidth.mean():.3f}, "
        f"coverage {g.fld_coverage_pct.mean():.1f}%")
wps = wpr.groupby("protein")[["mad_before", "mad_after", "loa_halfwidth_before", "loa_halfwidth_after"]].mean()
trs = trf.groupby("protein")[["mad_before", "mad_after", "loa_halfwidth_before", "loa_halfwidth_after"]].mean()
for nm, tb in (("within protein", wps), ("transfer from MS154", trs)):
    for p_, r in tb.iterrows():
        say(f"{nm}, {p_}: MAD {r.mad_before:.3f} -> {r.mad_after:.3f}, LoA {r.loa_halfwidth_before:.2f} -> {r.loa_halfwidth_after:.2f}")


# MS119 signal and replicate precision per protein (per-glycan files from lcms-compare fld)
import glob as _glob
for p_ in ("Trast", "Afli", "Cetux"):
    fs = sorted(_glob.glob(str(HERE / "data" / "ms119_per_glycan" / f"*{p_}*_fld_per_glycan.csv")))
    g = pd.concat([pd.read_csv(f).assign(run=Path(f).name) for f in fs])
    g = g[g.status == "matched"]
    tot = g.groupby("run")[["ms_area", "fld_area"]].sum()
    g["ms_p"] = 100 * g.ms_area / g.groupby("run").ms_area.transform("sum")
    g["fl_p"] = 100 * g.fld_area / g.groupby("run").fld_area.transform("sum")
    cv = g.groupby("glycan").agg(ms=("ms_p", lambda x: 100 * x.std(ddof=1) / x.mean()),
                                 fl=("fl_p", lambda x: 100 * x.std(ddof=1) / x.mean()))
    say(f"MS119 {p_}: {len(fs)} runs, MS/FLD summed area {(tot.ms_area / tot.fld_area).mean():.2f}, "
        f"replicate CV of relative abundance MS max {cv.ms.max():.1f}%, FLD max {cv.fl.max():.1f}%")

# Figure 1C: MS154 settings
for k, w in enumerate(sorted(sens.width.unique())):
    st = sens[sens.width == w].groupby("ppm").mean_abs_diff
    x = np.array(sorted(sens.ppm.unique())) * (1 + (k - 0.5) * 0.08)
    mean_, lo, hi = st.mean().values, st.min().values, st.max().values
    c1.errorbar(x, mean_, yerr=[mean_ - lo, hi - mean_], color=CAT[k], marker=MARK[k], ms=4.5, mec="white",
                mew=0.6, capsize=2, elinewidth=0.8, label=f"{w} points")
c1.set_xscale("log")
c1.set_xticks(sorted(sens.ppm.unique()))
c1.set_xticklabels([f"{p_:g}" for p_ in sorted(sens.ppm.unique())])
c1.minorticks_off()
c1.set_xlim(3.5, 28)
c1.set_ylim(0, 1.2)
c1.set_xlabel("Extraction window (± ppm)")
c1.set_ylabel("Mean |MS − FLD| (points)")
c1.yaxis.grid(True, color=GRID, lw=0.5)
c1.set_axisbelow(True)
c1.legend(title="Minimum peak width", title_fontsize=7, frameon=False, loc="upper right")
tag(c1, "C")
a1.set_title("2-AB fetuin, Q-TOF", fontsize=7, color=INK2)
b1.set_title("2-AB fetuin, Q-TOF", fontsize=7, color=INK2)
c1.set_title("RapiFluor-MS, Orbitrap", fontsize=7, color=INK2)
fig1.tight_layout(w_pad=1.6)
save(fig1, "Figure1_processing_settings")

# Figure 2: MS154 agreement
fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 2.6), gridspec_kw={"width_ratios": [1.05, 1.15, 1]})
y = np.arange(len(GLY))[::-1]
a.axvline(0, color=INK2, lw=0.6)
a.hlines(y, bias.lo, bias.hi, color=INK, lw=1.6)
a.plot(bias.bias, y, "o", color=INK, ms=4.5, mec="white", mew=0.6)
a.set_yticks(y)
a.set_yticklabels([GLAB[g] for g in GLY])
a.set_xlabel("MS − FLD (points), mean and 95% CI")
a.set_xlim(-1.2, 1.5)
a.xaxis.grid(True, color=GRID, lw=0.5)
a.set_axisbelow(True)
tag(a, "A", x=-0.38)
mean_ = (pr.ms_pct + pr.fld_pct) / 2
loa = 1.96 * dd.std(ddof=1)
for k, (day, g) in enumerate(pr.groupby("day")):
    b.plot((g.ms_pct + g.fld_pct) / 2, g.ms_pct - g.fld_pct, DAYM[day], color=CAT[k], ms=4, mec="white", mew=0.5,
           ls="", label=f"Day {day}")
for v, ls_ in ((0, "-"), (loa, "--"), (-loa, "--")):
    b.axhline(v, color=INK2, lw=0.7, ls=ls_)
b.text(0.95, loa, f"+{loa:.2f}", fontsize=6.3, color=INK2, va="bottom", ha="left")
b.text(0.95, -loa, f"−{loa:.2f}", fontsize=6.3, color=INK2, va="top", ha="left")
b.set_xscale("log")
b.set_xticks([1, 2, 5, 10, 20, 50])
b.set_xticklabels(["1", "2", "5", "10", "20", "50"])
b.minorticks_off()
b.set_xlim(0.9, 55)
b.set_ylim(-2, 2)
b.set_xlabel("Mean of MS and FLD (%)")
b.set_ylabel("MS − FLD (points)")
b.legend(frameon=False, loc="lower left", fontsize=6.3, handletextpad=0.1, ncol=3, columnspacing=0.8)
tag(b, "B")
c.axhline(1, color=INK2, lw=0.7, ls="--")
xg = np.arange(len(GLY))
for k, day in enumerate(sorted(rf_day.columns)):
    c.plot(xg + (k - 1) * 0.18, rf_day[day].values, DAYM[day], color=CAT[k], ms=4.2, mec="white", mew=0.5, ls="",
           label=f"Day {day}")
c.set_xticks(xg)
c.set_xticklabels([GLAB[g].replace(" (", "\n(") for g in GLY], fontsize=6.3)
c.set_ylabel("Response factor (MS% / FLD%)")
c.set_ylim(0.55, 1.2)
c.yaxis.grid(True, color=GRID, lw=0.5)
c.set_axisbelow(True)
c.legend(frameon=False, loc="lower right", fontsize=6.3, handletextpad=0.1)
tag(c, "C")
fig.tight_layout(w_pad=1.4)
save(fig, "Figure2_ms154_agreement")

# Figure 3: calibration tests
rows = [("Cross-day (MS154)", "Trastuzumab", cal.mad_before.mean(), cal.mad_after.mean())]
rows += [("Within protein (MS119)", PROT[p_], r.mad_before, r.mad_after) for p_, r in wps.loc[["Trast", "Afli", "Cetux"]].iterrows()]
rows += [("Transfer MS154 → MS119", PROT[p_], r.mad_before, r.mad_after) for p_, r in trs.loc[["Trast", "Afli", "Cetux"]].iterrows()]
fig, ax = plt.subplots(figsize=(3.4, 3.0))
ypos, labels, yy, prev = [], [], 0, None
for grp, prot, bef, aft in rows:
    if grp != prev:
        yy -= 0.45 if prev is not None else 0
        ax.text(0.052, yy, grp, fontsize=7, fontweight="bold", color=INK, va="center")
        yy -= 0.75
        prev = grp
    col = CAT[0] if aft < bef else CAT[1]
    ax.annotate("", xy=(aft, yy), xytext=(bef, yy), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2, shrinkA=3, shrinkB=3))
    ax.plot(bef, yy, "o", mfc="white", mec=col, mew=1.1, ms=5)
    ax.plot(aft, yy, "o", color=col, ms=5, mec="white", mew=0.5)
    ypos.append(yy)
    labels.append(prot)
    yy -= 0.75
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_xscale("log")
ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
ax.minorticks_off()
ax.set_xlim(0.05, 8)
ax.set_ylim(yy + 0.3, 0.5)
ax.set_xlabel("Mean |MS − FLD| (points)")
ax.xaxis.grid(True, color=GRID, lw=0.5)
ax.set_axisbelow(True)
ax.legend(handles=[Line2D([], [], ls="", marker="o", mfc="white", mec=INK2, ms=5, label="Before"),
                   Line2D([], [], ls="", marker="o", color=INK2, ms=5, label="After"),
                   Line2D([], [], color=CAT[0], lw=1.2, label="Decreased"),
                   Line2D([], [], color=CAT[1], lw=1.2, label="Increased")],
          frameon=False, fontsize=6.3, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=4, columnspacing=0.8,
          handletextpad=0.3)
fig.tight_layout()
save(fig, "Figure3_calibration")


# Figure S4: G1F isomer pairs (traces written by xic_pair.py)
PAIRS = [("pair_trast154_01.csv", "MS154 trastuzumab"), ("pair_trast119_01.csv", "MS119 trastuzumab"),
         ("pair_afli01.csv", "MS119 aflibercept")]
if all((HERE / "data" / "diag" / f).exists() for f, _ in PAIRS):
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.4), sharey=True)
    for k, (ax, (f, lab)) in enumerate(zip(axs, PAIRS)):
        tr = pd.read_csv(HERE / "data" / "diag" / f)
        sp = tr.split_min.iloc[0]
        txt = []
        for name, col, lw in (("FLD", INK2, 1.0), ("MS", CAT[0], 1.2)):
            q = tr[tr.trace == name]
            y = q.intensity.to_numpy()
            ax.plot(q.time_min, y / y.max(), color=col, lw=lw, label="Fluorescence" if name == "FLD" else "MS (G1F ions)")
            tq = q.time_min.to_numpy()
            l_ = __import__("scipy.integrate", fromlist=["trapezoid"]).trapezoid(y[tq <= sp], tq[tq <= sp])
            r_ = __import__("scipy.integrate", fromlist=["trapezoid"]).trapezoid(y[tq > sp], tq[tq > sp])
            txt.append(f"{name} {r_ / l_:.2f}")
            say(f"isomer pair {lab}: {name} late/early {r_ / l_:.3f}")
        ax.axvline(sp, color=INK2, lw=0.7, ls="--")
        ax.text(0.97, 0.95, "late/early\n" + "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=6.3, color=INK)
        ax.set_title(lab, fontsize=7, color=INK2)
        ax.set_xlabel("Retention time (min)")
        ax.set_ylim(-0.02, 1.08)
        tag(ax, "ABC"[k], x=-0.12 if k else -0.3)
    axs[0].set_ylabel("Normalized intensity")
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, frameon=False, fontsize=7, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(w_pad=1.0, rect=(0, 0, 1, 0.9))
    save(fig, "FigureS4_isomer_pairs")

# TOC graphic (3.25 x 1.75 in)
fig = plt.figure(figsize=(3.25, 1.75))
fig.text(0.03, 0.93, "GlyVAS: MS-only glycan profiles tested against fluorescence", fontsize=7, fontweight="bold", va="top")
t1 = fig.add_axes([0.12, 0.27, 0.36, 0.5])
t1.bar(range(len(GLY)), rf_all.values, color=CAT[0], width=0.65)
t1.axhline(1, color=INK2, lw=0.7, ls="--")
t1.set_xticks(range(len(GLY)))
t1.set_xticklabels([GLAB[g].replace(" (early)", "a").replace(" (late)", "b") for g in GLY], fontsize=5.5, rotation=40, ha="right", rotation_mode="anchor")
t1.set_ylim(0.5, 1.15)
t1.set_ylabel("MS / FLD", fontsize=6)
t1.tick_params(labelsize=5.5)
t1.set_title("Response factor", fontsize=6.5)
t2 = fig.add_axes([0.71, 0.27, 0.26, 0.5])
pairs = [("Same sequence", cal.mad_before.mean(), cal.mad_after.mean()),
         ("Next sequence", trs.loc["Trast", "mad_before"], trs.loc["Trast", "mad_after"])]
for i, (lab, bef, aft) in enumerate(pairs):
    col = CAT[0] if aft < bef else CAT[1]
    t2.annotate("", xy=(aft, -i), xytext=(bef, -i), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.1, shrinkA=2, shrinkB=2))
    t2.plot(bef, -i, "o", mfc="white", mec=col, ms=4)
    t2.plot(aft, -i, "o", color=col, ms=4, mec="white", mew=0.4)
t2.set_yticks([0, -1])
t2.set_yticklabels([p_[0] for p_ in pairs], fontsize=5.5)
t2.set_ylim(-1.6, 0.6)
t2.set_xlim(0, 0.8)
t2.tick_params(labelsize=5.5)
t2.set_xlabel("|MS − FLD| (points)", fontsize=6)
t2.set_title("After calibration", fontsize=6.5)
fig.savefig(OUT / "TOC_graphic.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
save(fig, "TOC_graphic")

(HERE / "numbers.txt").write_text("\n".join(lines) + "\n")
