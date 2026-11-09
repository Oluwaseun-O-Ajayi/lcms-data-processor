"""
Write a synthetic LC-MS run as mzML with known answers.

Two analytes elute as Gaussian peaks with known areas; a third ion with an m/z
close to analyte A (35 ppm away) co-elutes with it, so a narrow extraction
window must exclude it and a wide one includes it. Every MS1 spectrum also
carries random background ions. The file is written with pyOpenMS.

Usage:
    python examples/make_mzml_example.py [output.mzML]
"""

import sys
from pathlib import Path

import numpy as np

# Known answers. Areas are in intensity x minutes, the unit PeakDetector reports.
ANALYTES = {
    'A': {'mz': 1111.3870, 'rt_min': 6.0, 'sigma_min': 0.08, 'area': 2.0e5},
    'B': {'mz': 1257.4450, 'rt_min': 9.0, 'sigma_min': 0.10, 'area': 8.0e4},
}
# Interfering ion 35 ppm above analyte A, co-eluting with it.
INTERFERENT = {'mz': 1111.3870 * (1 + 35e-6), 'rt_min': 6.0, 'sigma_min': 0.08, 'area': 5.0e4}

RUN_MIN = 15.0
SCAN_INTERVAL_S = 0.5
SEED = 7


def elution(t_min, peak):
    return (peak['area'] / (peak['sigma_min'] * np.sqrt(2 * np.pi))
            * np.exp(-0.5 * ((t_min - peak['rt_min']) / peak['sigma_min']) ** 2))


def write_example(path):
    import pyopenms as oms

    rng = np.random.default_rng(SEED)
    exp = oms.MSExperiment()
    for rt_s in np.arange(0.0, RUN_MIN * 60, SCAN_INTERVAL_S):
        t = rt_s / 60.0
        mzs = list(rng.uniform(400, 2000, 40))
        ints = list(rng.exponential(300, 40))
        for peak in list(ANALYTES.values()) + [INTERFERENT]:
            mzs.append(peak['mz'])
            ints.append(elution(t, peak) + rng.normal(0, 50))
        order = np.argsort(mzs)
        spec = oms.MSSpectrum()
        spec.setRT(float(rt_s))
        spec.setMSLevel(1)
        spec.set_peaks((np.asarray(mzs)[order], np.clip(np.asarray(ints)[order], 0, None)))
        exp.addSpectrum(spec)
    oms.MzMLFile().store(str(path), exp)
    return Path(path)


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'data' / 'example_run.mzML'
    write_example(out)
    print(f"Wrote {out}")
