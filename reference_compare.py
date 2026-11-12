"""
Compare MS-based glycan quantitation with independent references.

Two references are supported:

1. Fluorescence (FLD). For 2-AB and other fluorescent labels, every glycan carries one
   label, so fluorescence peak area is proportional to molar amount. The FLD trace is
   read from the same mzML file (MSConvert exports detector traces as chromatograms).
   The FLD-to-MS transfer delay is estimated by cross-correlating the FLD trace with the
   MS total-ion chromatogram unless it is given. Each glycan is matched to the FLD peak
   whose apex lies closest to its MS apex minus the delay. Glycans that fall on the same
   FLD peak (co-elution) cannot be separated by fluorescence and are reported but left
   out of the comparison.

2. Skyline. `skyline-targets` writes a small-molecule transition list from the same
   target m/z values; after integrating in Skyline and exporting a report, `skyline`
   compares Skyline's areas with this package's areas.

Agreement is reported as the mean absolute difference in relative abundance
(percentage points), Bland-Altman bias with 95% limits of agreement, and Pearson
correlation; for Skyline also the Spearman correlation and median ratio of areas.

Usage:
    lcms-compare fld --mzml run.mzML --targets results/run_targets.csv --fld FLD \\
                     --ppm 10 --width 5 --baseline-window 0.1 --out results/run_fld
    lcms-compare skyline-targets --targets results/run_target_list.csv --out skyline_targets.csv
    lcms-compare skyline --targets results/run_targets.csv --skyline skyline_report.csv \\
                     --ppm 10 --width 5 --baseline-window 0.1 --out results/run_skyline
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from lcms_data_processor import MzMLReader, PeakDetector


# --- agreement statistics ---------------------------------------------------

def agreement(a, b):
    """Agreement between two sets of relative abundances (percent)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = a - b
    out = {'n': len(a), 'mean_abs_diff_pct_points': float(np.mean(np.abs(d))),
           'bias_pct_points': float(np.mean(d)),
           'loa_low': float(np.mean(d) - 1.96 * np.std(d, ddof=1)) if len(d) > 1 else np.nan,
           'loa_high': float(np.mean(d) + 1.96 * np.std(d, ddof=1)) if len(d) > 1 else np.nan,
           'pearson_r': float(stats.pearsonr(a, b)[0]) if len(a) > 2 else np.nan}
    return out


def pick_setting(df, ppm, width, baseline_window):
    sel = df[np.isclose(df['ppm'], ppm) & (df['width'] == width) & np.isclose(df['baseline_window'], baseline_window)]
    if sel.empty:
        raise SystemExit(f"No rows for ppm={ppm}, width={width}, baseline_window={baseline_window}")
    return sel


# --- fluorescence -----------------------------------------------------------

def estimate_delay(t1, y1, t2, y2, max_lag=1.0, step=None):
    """
    Lag (minutes) that best aligns trace 2 to trace 1: y1(t) ~ y2(t - lag).
    A positive lag means trace 1 is later (for example MS after FLD).
    """
    step = step or float(np.median(np.diff(t1)))
    grid = np.arange(max(t1.min(), t2.min()) + max_lag, min(t1.max(), t2.max()) - max_lag, step)
    a = np.interp(grid, t1, y1)
    a = (a - a.mean()) / (a.std() or 1)
    lags = np.arange(-max_lag, max_lag + step / 2, step)
    scores = []
    for lag in lags:
        b = np.interp(grid - lag, t2, y2)
        b = (b - b.mean()) / (b.std() or 1)
        scores.append(float(np.mean(a * b)))
    return float(lags[int(np.argmax(scores))]), float(max(scores))


def ms_glycan_table(targets, setting):
    """Per glycan: MS area summed over charge states and apex time of the strongest charge state."""
    sel = setting[setting['peak_area'] > 0]
    area = sel.groupby('glycan')['peak_area'].sum()
    apex = sel.loc[sel.groupby('glycan')['peak_area'].idxmax()].set_index('glycan')['retention_time']
    return pd.DataFrame({'ms_area': area, 'ms_rt': apex}).reset_index()


def ms_apexes(setting, merge=0.05):
    """Apex times per glycan across charge states (from 'peak_rts' when present), merged within `merge` min."""
    out = {}
    sel = setting[setting['peak_area'] > 0]
    for g, grp in sel.groupby('glycan'):
        if 'peak_rts' in grp and grp['peak_rts'].notna().any():
            rts = [float(x) for v in grp['peak_rts'].dropna().astype(str) for x in v.split(';') if x]
        else:
            rts = grp['retention_time'].dropna().astype(float).tolist()
        merged = []
        for t in sorted(rts):
            if not merged or t - merged[-1] > merge:
                merged.append(t)
        out[g] = merged
    return out


def compare_fld(reader, targets_long, fld_match, ppm, width, baseline_window,
                delay=None, tolerance=0.1, max_lag=1.0, fld_width=5, fld_baseline_window=0.1,
                fld_prominence=0.02, margin=2.0, return_peaks=False):
    setting = pick_setting(targets_long, ppm, width, baseline_window)
    ms = ms_glycan_table(targets_long, setting)
    fld = reader.chromatogram(fld_match)

    if delay is None:
        tic = reader.tic()
        delay, score = estimate_delay(tic['time'].to_numpy(), tic['intensity'].to_numpy(),
                                      fld['time'].to_numpy(), fld['intensity'].to_numpy(), max_lag=max_lag)
    else:
        score = np.nan

    # Detect fluorescence peaks only where the target glycans elute (plus a margin), so the
    # free-label peak and other large early signals do not set the detection threshold.
    lo = ms['ms_rt'].min() - delay - margin
    hi = ms['ms_rt'].max() - delay + margin
    region = fld[(fld['time'] >= lo) & (fld['time'] <= hi)].reset_index(drop=True)
    det = PeakDetector(region)
    det.detect_peaks(prominence=fld_prominence, width=fld_width, baseline_window=fld_baseline_window)
    det.integrate_peaks()
    peaks = det.peaks.reset_index(drop=True)

    # Every MS apex of a glycan (several when isomer peaks were summed) is matched to the
    # nearest FLD peak within the tolerance; the glycan's FLD area is the sum of its FLD peaks.
    apexes = ms_apexes(setting)
    ms['expected_fld_rt'] = ms['ms_rt'] - delay
    claims, fld_ids, fld_rts, fld_area = {}, [], [], []
    for _, row in ms.iterrows():
        ids = set()
        for rt in apexes.get(row['glycan'], [row['ms_rt']]):
            if peaks.empty:
                break
            dist = (peaks['retention_time'] - (rt - delay)).abs()
            j = dist.idxmin()
            if dist[j] <= tolerance:
                ids.add(int(peaks.loc[j, 'peak_id']))
        for k in ids:
            claims[k] = claims.get(k, 0) + 1
        sel = peaks[peaks['peak_id'].isin(ids)]
        fld_ids.append(';'.join(str(k) for k in sorted(ids)))
        fld_rts.append(';'.join(f"{t:.3f}" for t in sorted(sel['retention_time'])))
        fld_area.append(sel['peak_area'].sum() if ids else np.nan)
    ms['fld_peak_ids'], ms['fld_rts'], ms['fld_area'] = fld_ids, fld_rts, fld_area
    shared = [any(claims[int(k)] > 1 for k in ids.split(';')) if ids else False for ids in fld_ids]
    ms['status'] = np.where(ms['fld_peak_ids'] == '', 'no FLD peak',
                            np.where(shared, 'shared FLD peak', 'matched'))
    ok = ms['status'] == 'matched'
    ms['ms_pct'] = np.nan
    ms['fld_pct'] = np.nan
    if ok.any():
        ms.loc[ok, 'ms_pct'] = 100 * ms.loc[ok, 'ms_area'] / ms.loc[ok, 'ms_area'].sum()
        ms.loc[ok, 'fld_pct'] = 100 * ms.loc[ok, 'fld_area'] / ms.loc[ok, 'fld_area'].sum()
    summary = agreement(ms.loc[ok, 'ms_pct'], ms.loc[ok, 'fld_pct']) if ok.sum() else {'n': 0}
    summary.update({'delay_min': delay, 'delay_correlation': score, 'ppm': ppm, 'width': width,
                    'baseline_window': baseline_window, 'fld_region_min': lo, 'fld_region_max': hi,
                    'fld_peaks_detected': len(peaks)})
    if return_peaks:
        return ms, summary, peaks
    return ms, summary


# --- Skyline ----------------------------------------------------------------

def skyline_transition_list(target_list):
    """Precursor-only small-molecule transition list for Skyline (MS1 full-scan filtering)."""
    t = target_list.copy()
    if 'charge' not in t:
        raise SystemExit("target list needs a 'charge' column (use the _target_list.csv from lcms-batch --glycans)")
    out = pd.DataFrame({'Molecule List Name': 'glycans', 'Molecule Name': t['glycan'],
                        'Precursor m/z': t['mz'].round(5), 'Precursor Charge': t['charge'].astype(int)})
    for col, src in (('Explicit Retention Time', 'rt_min'),):
        if src in t and t[src].notna().any():
            out[col] = (t['rt_min'] + t['rt_max']) / 2
            out['Explicit Retention Time Window'] = t['rt_max'] - t['rt_min']
    return out


AREA_COLUMNS = ('Total Area MS1', 'Total Area', 'Area')


def compare_skyline(targets_long, report, ppm, width, baseline_window, area_col=None):
    setting = pick_setting(targets_long, ppm, width, baseline_window).copy()
    if 'charge' not in setting or setting['charge'].isna().all():
        setting['charge'] = setting['target'].str.extract(r'_(\d+)\+$')[0].astype(float)
    area_col = area_col or next((c for c in AREA_COLUMNS if c in report), None)
    need = {'Molecule Name', 'Precursor Charge'}
    if area_col is None or not need <= set(report.columns):
        raise SystemExit(f"Skyline report needs columns {sorted(need)} and one of {AREA_COLUMNS}")
    sky = report.rename(columns={'Molecule Name': 'glycan', 'Precursor Charge': 'charge', area_col: 'skyline_area'})
    sky['skyline_area'] = pd.to_numeric(sky['skyline_area'], errors='coerce')
    sky = sky.groupby(['glycan', 'charge'], as_index=False)['skyline_area'].sum()
    m = setting.merge(sky, on=['glycan', 'charge'], how='inner')
    both = m[(m['peak_area'] > 0) & (m['skyline_area'] > 0)]

    per_glycan = both.groupby('glycan')[['peak_area', 'skyline_area']].sum()
    per_glycan['this_pct'] = 100 * per_glycan['peak_area'] / per_glycan['peak_area'].sum()
    per_glycan['skyline_pct'] = 100 * per_glycan['skyline_area'] / per_glycan['skyline_area'].sum()
    summary = agreement(per_glycan['this_pct'], per_glycan['skyline_pct'])
    if len(both) > 2:
        summary['spearman_rho_areas'] = float(stats.spearmanr(both['peak_area'], both['skyline_area'])[0])
        summary['pearson_r_log_areas'] = float(stats.pearsonr(np.log10(both['peak_area']),
                                                              np.log10(both['skyline_area']))[0])
    summary['median_area_ratio'] = float(np.median(both['peak_area'] / both['skyline_area'])) if len(both) else np.nan
    summary['n_precursors'] = len(both)
    return m, per_glycan.reset_index(), summary


# --- command line -----------------------------------------------------------

def _setting_args(p):
    p.add_argument('--ppm', type=float, default=10.0)
    p.add_argument('--width', type=int, default=5)
    p.add_argument('--baseline-window', type=float, default=0.1)


def build_parser():
    p = argparse.ArgumentParser(prog='lcms-compare', description=__doc__.split('\n\n')[0])
    sub = p.add_subparsers(dest='cmd', required=True)

    f = sub.add_parser('fld', help='compare MS relative abundances with fluorescence')
    f.add_argument('--mzml', required=True)
    f.add_argument('--targets', required=True, help='<run>_targets.csv from lcms-batch')
    f.add_argument('--fld', default='FLD', help='index or native-ID substring of the FLD chromatogram')
    f.add_argument('--delay', type=float, help='FLD-to-MS delay in minutes (default: estimate)')
    f.add_argument('--tolerance', type=float, default=0.1, help='max RT difference for matching (min)')
    f.add_argument('--max-lag', type=float, default=1.0, help='delay search range (min)')
    f.add_argument('--fld-prominence', type=float, default=0.02,
                   help='FLD peak prominence as a fraction of the largest peak in the glycan region')
    f.add_argument('--fld-width', type=int, default=5, help='minimum FLD peak width (points)')
    f.add_argument('--margin', type=float, default=2.0,
                   help='FLD region extends this many minutes beyond the first and last glycan')
    f.add_argument('--list', action='store_true', help='only list the chromatograms in the file')
    f.add_argument('--out', default='results/fld')
    _setting_args(f)

    st = sub.add_parser('skyline-targets', help='write a Skyline transition list')
    st.add_argument('--targets', required=True, help='<run>_target_list.csv from lcms-batch')
    st.add_argument('--out', default='skyline_transition_list.csv')

    s = sub.add_parser('skyline', help='compare with a Skyline report')
    s.add_argument('--targets', required=True, help='<run>_targets.csv from lcms-batch')
    s.add_argument('--skyline', required=True, help='Skyline report CSV')
    s.add_argument('--area-column', help=f'area column in the report (default: first of {AREA_COLUMNS})')
    s.add_argument('--out', default='results/skyline')
    _setting_args(s)
    return p


def _write_summary(summary, path):
    pd.DataFrame([summary]).to_csv(path, index=False)
    for k, v in summary.items():
        print(f"  {k}: {v:.4g}" if isinstance(v, float) else f"  {k}: {v}")


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.cmd == 'skyline-targets':
        skyline_transition_list(pd.read_csv(a.targets)).to_csv(a.out, index=False)
        print(f"Wrote {a.out}")
        return
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if a.cmd == 'fld':
        reader = MzMLReader(a.mzml)
        if a.list:
            print(reader.list_chromatograms().to_string(index=False))
            return
        fld_match = int(a.fld) if a.fld.isdigit() else a.fld
        table, summary, peaks = compare_fld(reader, pd.read_csv(a.targets), fld_match, a.ppm, a.width,
                                            a.baseline_window, a.delay, a.tolerance, a.max_lag,
                                            fld_width=a.fld_width, fld_prominence=a.fld_prominence,
                                            margin=a.margin, return_peaks=True)
        table.to_csv(f"{out}_per_glycan.csv", index=False)
        peaks.to_csv(f"{out}_fld_peaks.csv", index=False)
        print(table[['glycan', 'ms_rt', 'fld_rts', 'status', 'ms_pct', 'fld_pct']].to_string(index=False))
        _write_summary(summary, f"{out}_summary.csv")
    else:
        merged, per_glycan, summary = compare_skyline(pd.read_csv(a.targets), pd.read_csv(a.skyline),
                                                      a.ppm, a.width, a.baseline_window, a.area_column)
        merged.to_csv(f"{out}_per_precursor.csv", index=False)
        per_glycan.to_csv(f"{out}_per_glycan.csv", index=False)
        print(per_glycan.to_string(index=False))
        _write_summary(summary, f"{out}_summary.csv")


if __name__ == '__main__':
    sys.exit(main())
