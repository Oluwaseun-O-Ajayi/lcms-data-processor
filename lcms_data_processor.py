"""
LC-MS Data Processing Pipeline
===============================

A Python toolkit for processing chromatograms exported from liquid
chromatography-mass spectrometry (LC-MS) instruments as time-intensity CSV
files: baseline estimation, peak detection and integration, signal-to-noise,
peak-shape and resolution checks, and linear calibration.

Features:
- Chromatogram data import and visualization
- Automated peak detection and integration
- Calibration curve generation and quantification
- Quality control metrics (S/N ratio, peak shape, resolution)
- Batch sample processing
- Publication-quality reports

Author: Oluwaseun O. Ajayi
Institution: University of Georgia
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal, integrate
from pathlib import Path


class ChromatogramParser:
    """
    Parse and extract information from chromatogram data files.
    """
    
    def __init__(self, data_file):
        """
        Initialize chromatogram parser.
        
        Args:
            data_file: Path to CSV file with 'time' and 'intensity' columns
        """
        self.data_file = Path(data_file)
        self.data = None
        self.metadata = {}
        
        if self.data_file.exists():
            self._load_data()
        else:
            raise FileNotFoundError(f"Data file not found: {data_file}")
    
    def _load_data(self):
        """Load chromatogram data from CSV file."""
        try:
            self.data = pd.read_csv(self.data_file)
            
            # Validate required columns
            required_cols = ['time', 'intensity']
            if not all(col in self.data.columns for col in required_cols):
                raise ValueError(f"CSV must contain columns: {required_cols}")
            
            # Sort by time
            self.data = self.data.sort_values('time').reset_index(drop=True)
            
            # Extract metadata from filename if present
            self.metadata['filename'] = self.data_file.name
            self.metadata['num_points'] = len(self.data)
            self.metadata['time_range'] = (self.data['time'].min(), 
                                          self.data['time'].max())
            self.metadata['max_intensity'] = self.data['intensity'].max()
            
        except Exception as e:
            raise Exception(f"Error loading chromatogram data: {e}")
    
    def get_chromatogram(self):
        """Get chromatogram data as DataFrame."""
        return self.data.copy()
    
    def plot_chromatogram(self, save_path=None):
        """Plot raw chromatogram."""
        fig, ax = plt.subplots(figsize=(12, 6))
        
        ax.plot(self.data['time'], self.data['intensity'], 
               linewidth=1.5, color='steelblue', alpha=0.8)
        
        ax.set_xlabel('Retention Time (min)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Intensity (AU)', fontsize=12, fontweight='bold')
        ax.set_title(f'Chromatogram: {self.metadata["filename"]}', 
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add metadata text
        info_text = (f"Points: {self.metadata['num_points']}\n"
                    f"Time: {self.metadata['time_range'][0]:.2f}-"
                    f"{self.metadata['time_range'][1]:.2f} min\n"
                    f"Max Intensity: {self.metadata['max_intensity']:.2e}")
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✅ Chromatogram plot saved to: {save_path}")
        
        plt.close(fig)


class PeakDetector:
    """
    Detect and integrate peaks in chromatograms.
    """
    
    def __init__(self, chromatogram_data):
        """
        Initialize peak detector.
        
        Args:
            chromatogram_data: DataFrame with 'time' and 'intensity' columns
        """
        self.data = chromatogram_data.copy()
        self.peaks = None
        self.baseline = None
    
    def detect_peaks(self, prominence=0.1, width=5, height=None, baseline_window=0.1,
                     limit_halfwidths=3.4):
        """
        Detect peaks with scipy.signal.find_peaks on the baseline-corrected trace.

        Args:
            prominence: minimum prominence as a fraction of the largest corrected intensity
                (never below three noise standard deviations)
            width: minimum peak width in data points
            height: minimum corrected height (absolute); default is the larger of 5% of the
                largest corrected intensity and three noise standard deviations
            baseline_window: baseline window as a fraction of the number of points

        Each peak records its apex index, retention time, height, full width at half
        maximum (FWHM, in time units) and the integration limits. Each limit lies
        `limit_halfwidths` half-widths at half maximum from the apex on its own side
        (3.4 half-widths correspond to four standard deviations of a Gaussian peak), so
        tailing is followed on the trailing side. Limits never cross the valley between
        two neighbouring detected peaks.
        """
        self.baseline = self._estimate_baseline(baseline_window)
        corrected = (self.data['intensity'] - self.baseline).to_numpy()
        self.noise_sd = self._estimate_noise(corrected)

        # Detection and widths use a Savitzky-Golay smoothed trace (window 2*width+1
        # points, quadratic) so that noise spikes do not split peaks or shrink widths;
        # heights and areas are taken from the unsmoothed corrected trace.
        window = 2 * int(width) + 1
        smooth = signal.savgol_filter(corrected, window, 2) if len(corrected) > window else corrected
        if height is None:
            height = max(0.05 * smooth.max(), 3.0 * self.noise_sd)
        idx, props = signal.find_peaks(smooth, prominence=max(prominence * smooth.max(), 3.0 * self.noise_sd),
                                       width=width, height=height)
        half = signal.peak_widths(smooth, idx, rel_height=0.5)
        t = self.data['time'].to_numpy()
        dt = float(np.median(np.diff(t))) if len(t) > 1 else 1.0

        rows = []
        n = len(smooth)
        for i, apex in enumerate(idx):
            # Limits never cross the valley between neighbouring detected peaks.
            lo = int(idx[i - 1] + np.argmin(smooth[idx[i - 1]:apex])) if i > 0 else 0
            hi = int(apex + np.argmin(smooth[apex:idx[i + 1] + 1])) if i < len(idx) - 1 else n - 1
            hw_left, hw_right = apex - half[2][i], half[3][i] - apex
            left = max(lo, int(np.floor(apex - limit_halfwidths * hw_left)))
            right = min(hi, int(np.ceil(apex + limit_halfwidths * hw_right)))
            rows.append({
                'peak_id': i + 1,
                'apex_index': int(apex),
                'retention_time': float(t[apex]),
                'peak_height': float(smooth[apex]),
                'prominence': float(props['prominences'][i]),
                'width': float(half[0][i]),            # FWHM in data points
                'fwhm': float(half[0][i] * dt),        # FWHM in time units
                'left_base': left,
                'right_base': right,
            })
        self.peaks = pd.DataFrame(rows, columns=['peak_id', 'apex_index', 'retention_time', 'peak_height',
                                                 'prominence', 'width', 'fwhm', 'left_base', 'right_base'])

        # Second pass: re-estimate the baseline from the points outside the peak windows
        # (rolling median, linear interpolation across the windows). This keeps broad
        # peaks from lifting the baseline estimate underneath themselves.
        if len(self.peaks):
            y = self.data['intensity'].to_numpy(dtype=float)
            mask = np.ones(n, dtype=bool)
            for p in self.peaks.itertuples():
                mask[int(p.left_base):int(p.right_base) + 1] = False
            if mask.sum() >= 2:
                outside = pd.Series(np.where(mask, y, np.nan))
                med = outside.rolling(window, center=True, min_periods=1).median().to_numpy()
                ok = mask & np.isfinite(med)
                base2 = np.interp(np.arange(n), np.flatnonzero(ok), med[ok])
                self.baseline = pd.Series(base2, index=self.data.index)
                corrected = y - base2
                smooth = signal.savgol_filter(corrected, window, 2) if n > window else corrected
                self.peaks['peak_height'] = [float(smooth[a]) for a in self.peaks['apex_index']]
        return self.peaks

    def _estimate_baseline(self, window_fraction=0.1):
        """
        Baseline by morphological opening (rolling minimum followed by rolling maximum)
        over a window of window_fraction of the points, smoothed with a rolling mean.

        An opening follows the lower envelope of the noise, so it lies below the true
        baseline by roughly two to three noise standard deviations. The offset is
        corrected by adding the median of (intensity - opening) over points that are
        within the noise band, which removes the positive area bias that a plain
        rolling minimum produces.
        """
        y = self.data['intensity']
        window = max(int(len(y) * window_fraction), 3)
        opening = (y.rolling(window, center=True, min_periods=1).min()
                    .rolling(window, center=True, min_periods=1).max()
                    .rolling(window, center=True, min_periods=1).mean())
        residual = (y - opening).to_numpy()
        sd = self._estimate_noise(residual)
        quiet = residual < np.median(residual) + 3 * sd
        offset = float(np.median(residual[quiet])) if quiet.any() else 0.0
        return opening + offset

    @staticmethod
    def _estimate_noise(values):
        """Robust noise standard deviation from first differences:
        1.4826 * MAD(diff) / sqrt(2). Smooth peaks contribute little to first differences."""
        d = np.diff(np.asarray(values, dtype=float))
        if len(d) == 0:
            return 0.0
        mad = np.median(np.abs(d - np.median(d)))
        return float(1.4826 * mad / np.sqrt(2.0))

    def integrate_peaks(self):
        """Integrate each peak over its integration limits with the trapezoidal rule
        (area in intensity x time units)."""
        if self.peaks is None:
            raise ValueError("No peaks detected. Run detect_peaks() first.")
        corrected = (self.data['intensity'] - self.baseline).to_numpy()
        t = self.data['time'].to_numpy()
        self.peaks['peak_area'] = [
            float(integrate.trapezoid(corrected[int(p.left_base):int(p.right_base) + 1],
                                      t[int(p.left_base):int(p.right_base) + 1]))
            for p in self.peaks.itertuples()
        ]
        return self.peaks

    def calculate_snr(self):
        """Signal-to-noise ratio as corrected peak height divided by the noise standard
        deviation. Note that pharmacopeial definitions (for example 2H/h with peak-to-peak
        noise) give different numbers."""
        if self.peaks is None:
            raise ValueError("No peaks detected. Run detect_peaks() first.")
        self.peaks['snr'] = self.peaks['peak_height'] / max(self.noise_sd, 1e-12)
        return self.peaks

    def plot_detected_peaks(self, save_path=None):
        """Visualize detected peaks on chromatogram."""
        if self.peaks is None:
            raise ValueError("No peaks detected. Run detect_peaks() first.")
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
        
        # Plot 1: Chromatogram with peaks
        corrected_intensity = self.data['intensity'] - self.baseline
        
        ax1.plot(self.data['time'], corrected_intensity, 
                linewidth=1.5, color='steelblue', alpha=0.7, 
                label='Baseline-corrected')
        ax1.plot(self.data['time'], self.baseline, 
                linewidth=1, color='gray', linestyle='--', 
                label='Baseline')
        
        # Mark detected peaks
        for _, peak in self.peaks.iterrows():
            ax1.plot(peak['retention_time'], peak['peak_height'],
                    'ro', markersize=10)
            ax1.axvline(peak['retention_time'], color='red', 
                       linestyle=':', alpha=0.3)
            ax1.text(peak['retention_time'], peak['peak_height'],
                    f"  {peak['peak_id']}", fontsize=10, fontweight='bold')
        
        ax1.set_xlabel('Retention Time (min)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Intensity (AU)', fontsize=11, fontweight='bold')
        ax1.set_title('Detected Peaks', fontsize=13, fontweight='bold')
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Peak metrics
        x_pos = np.arange(len(self.peaks))

        ax2_twin = ax2.twinx()

        bars1 = ax2.bar(x_pos - 0.2, self.peaks['peak_area'], 
               width=0.4, label='Peak Area', 
               color='skyblue', edgecolor='black')
        bars2 = ax2_twin.bar(x_pos + 0.2, self.peaks['snr'], 
                    width=0.4, label='S/N Ratio',
                    color='coral', edgecolor='black')

        ax2.set_xlabel('Peak ID', fontsize=11, fontweight='bold', labelpad=10)
        ax2.set_ylabel('Peak Area', fontsize=11, fontweight='bold', color='skyblue', labelpad=10)
        ax2_twin.set_ylabel('S/N Ratio', fontsize=11, fontweight='bold', color='coral', labelpad=10)
        ax2.set_title('Peak Metrics', fontsize=13, fontweight='bold', pad=15)
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(self.peaks['peak_id'], fontsize=10)
        ax2.grid(True, alpha=0.3, axis='y')

        # Color the y-axis labels to match the bars
        ax2.tick_params(axis='y', labelcolor='skyblue', labelsize=10)
        ax2_twin.tick_params(axis='y', labelcolor='coral', labelsize=10)
        
        # Combine legends
        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2_twin.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)
        
        plt.tight_layout(pad=2.0)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✅ Peak detection plot saved to: {save_path}")
        
        plt.close(fig)
        
        # Print peak summary
        print("\n=== Peak Detection Summary ===")
        print(f"Total peaks detected: {len(self.peaks)}")
        print("\nPeak Details:")
        print(self.peaks[['peak_id', 'retention_time', 'peak_area', 'snr']].to_string(index=False))


class CalibrationCurve:
    """
    Generate calibration curves and quantify analytes.
    """
    
    def __init__(self):
        """Initialize calibration curve."""
        self.standards = None
        self.curve_params = None
        self.r_squared = None
    
    def fit_calibration(self, concentrations, peak_areas):
        """
        Fit calibration curve using linear regression.
        
        Args:
            concentrations: Array of standard concentrations
            peak_areas: Array of corresponding peak areas
        """
        concentrations = np.array(concentrations)
        peak_areas = np.array(peak_areas)
        
        # Store standards
        self.standards = pd.DataFrame({
            'concentration': concentrations,
            'peak_area': peak_areas
        })
        
        if len(concentrations) < 3:
            raise ValueError("At least three calibration standards are required")
        # Fit linear model: area = slope * concentration + intercept
        self.curve_params = np.polyfit(concentrations, peak_areas, 1)
        
        # Calculate R²
        predicted_areas = np.polyval(self.curve_params, concentrations)
        ss_res = np.sum((peak_areas - predicted_areas) ** 2)
        ss_tot = np.sum((peak_areas - np.mean(peak_areas)) ** 2)
        self.r_squared = 1 - (ss_res / ss_tot)
        
        return self.curve_params, self.r_squared

    @property
    def slope(self):
        return None if self.curve_params is None else float(self.curve_params[0])

    @property
    def intercept(self):
        return None if self.curve_params is None else float(self.curve_params[1])
    
    def quantify_samples(self, peak_areas):
        """
        Quantify unknown samples using calibration curve.
        
        Args:
            peak_areas: Array or single value of peak areas
        """
        if self.curve_params is None:
            raise ValueError("No calibration curve fitted. Run fit_calibration() first.")
        
        # Convert to array if single value
        peak_areas = np.atleast_1d(peak_areas)
        
        # Calculate concentrations: c = (y - b) / m
        concentrations = (peak_areas - self.curve_params[1]) / self.curve_params[0]
        
        return concentrations
    
    def plot_calibration_curve(self, save_path=None):
        """Plot calibration curve with standards."""
        if self.standards is None:
            raise ValueError("No calibration data. Run fit_calibration() first.")
        
        fig, ax = plt.subplots(figsize=(10, 7))
        
        # Plot standards
        ax.scatter(self.standards['concentration'], 
                  self.standards['peak_area'],
                  s=150, color='darkblue', edgecolors='black', 
                  linewidth=2, label='Standards', zorder=3)
        
        # Plot fitted line
        conc_range = np.linspace(0, self.standards['concentration'].max() * 1.1, 100)
        fitted_areas = np.polyval(self.curve_params, conc_range)
        ax.plot(conc_range, fitted_areas, 'r--', linewidth=2, 
               label='Calibration Curve', zorder=2)
        
        # Add equation and R²
        slope, intercept = self.curve_params
        equation = f'y = {slope:.2e}x + {intercept:.2e}'
        r_sq_text = f'R² = {self.r_squared:.4f}'
        
        ax.text(0.05, 0.95, f'{equation}\n{r_sq_text}',
               transform=ax.transAxes, fontsize=12,
               verticalalignment='top', fontweight='bold',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        ax.set_xlabel('Concentration', fontsize=12, fontweight='bold')
        ax.set_ylabel('Peak Area', fontsize=12, fontweight='bold')
        ax.set_title('Calibration Curve', fontsize=14, fontweight='bold')
        ax.legend(loc='lower right', fontsize=11)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✅ Calibration curve plot saved to: {save_path}")
        
        plt.close(fig)
        
        # Print calibration info
        print("\n=== Calibration Curve Summary ===")
        print(f"Slope: {slope:.4e}")
        print(f"Intercept: {intercept:.4e}")
        print(f"R²: {self.r_squared:.4f}")
        
        if self.r_squared < 0.95:
            print("⚠️ Warning: R² < 0.95 - calibration may be poor")
        else:
            print("✅ Good calibration (R² ≥ 0.95)")


class QualityControl:
    """
    Perform quality control checks on LC-MS data.
    """
    
    def __init__(self):
        """Initialize QC checker."""
        self.qc_results = {}
    
    def check_peak_shape(self, peak_data, chromatogram_data, baseline=None, smooth_window=11):
        """
        Asymmetry factor As = b / a measured at 10% of peak height, where a and b are the
        distances from the apex to the leading and trailing edges. As = 1 for a symmetric
        peak and As > 1 for tailing.

        Args:
            peak_data: peaks from PeakDetector (needs 'apex_index' and 'peak_id')
            chromatogram_data: DataFrame with 'time' and 'intensity'
            baseline: baseline to subtract (PeakDetector.baseline); if None, the minimum
                intensity is used
            smooth_window: Savitzky-Golay window (points) applied before measuring widths
        """
        y = chromatogram_data['intensity'].to_numpy(dtype=float)
        y = y - (np.asarray(baseline, dtype=float) if baseline is not None else y.min())
        if len(y) > smooth_window:
            y = signal.savgol_filter(y, smooth_window, 2)
        apexes = peak_data['apex_index'].to_numpy(dtype=int)
        rows = []
        if len(apexes):
            _, _, left_ips, right_ips = signal.peak_widths(y, apexes, rel_height=0.9)
            for pid, apex, l, r in zip(peak_data['peak_id'], apexes, left_ips, right_ips):
                a, b = apex - l, r - apex
                asym = b / a if a > 0 else np.nan
                rows.append({'peak_id': pid, 'asymmetry_factor': asym,
                             'quality': 'Good' if 0.9 <= asym <= 1.2 else
                                        'Acceptable' if 0.8 <= asym <= 1.5 else 'Poor'})
        self.qc_results['peak_shape'] = pd.DataFrame(rows, columns=['peak_id', 'asymmetry_factor', 'quality'])
        return self.qc_results['peak_shape']

    def check_resolution(self, peak_data):
        """
        Resolution of adjacent peaks from half-height widths:
        Rs = 1.18 (t2 - t1) / (w0.5,1 + w0.5,2). Rs >= 1.5 indicates baseline separation.

        Args:
            peak_data: peaks from PeakDetector (needs 'retention_time' and 'fwhm')
        """
        if len(peak_data) < 2:
            self.qc_results['resolution'] = pd.DataFrame(columns=['peak_pair', 'resolution', 'quality'])
            return self.qc_results['resolution']
        ordered = peak_data.sort_values('retention_time').reset_index(drop=True)
        rows = []
        for i in range(len(ordered) - 1):
            p1, p2 = ordered.iloc[i], ordered.iloc[i + 1]
            denom = p1['fwhm'] + p2['fwhm']
            rs = 1.18 * (p2['retention_time'] - p1['retention_time']) / denom if denom > 0 else np.nan
            rows.append({'peak_pair': f"{int(p1['peak_id'])}-{int(p2['peak_id'])}", 'resolution': rs,
                         'quality': 'Baseline' if rs >= 1.5 else 'Partial' if rs >= 1.0 else 'Poor'})
        self.qc_results['resolution'] = pd.DataFrame(rows)
        return self.qc_results['resolution']

    def generate_qc_report(self, save_path=None):
        """Generate comprehensive QC report."""
        print("\n" + "="*70)
        print("   QUALITY CONTROL REPORT")
        print("="*70)
        
        if 'peak_shape' in self.qc_results:
            print("\n📊 Peak Shape Analysis:")
            print(self.qc_results['peak_shape'].to_string(index=False))
            
            good_peaks = (self.qc_results['peak_shape']['quality'] == 'Good').sum()
            total_peaks = len(self.qc_results['peak_shape'])
            print(f"\nPeak Shape Summary: {good_peaks}/{total_peaks} peaks with good symmetry")
        
        if 'resolution' in self.qc_results:
            print("\n🔍 Peak Resolution Analysis:")
            print(self.qc_results['resolution'].to_string(index=False))
            
            baseline_sep = (self.qc_results['resolution']['quality'] == 'Baseline').sum()
            total_pairs = len(self.qc_results['resolution'])
            print(f"\nResolution Summary: {baseline_sep}/{total_pairs} peak pairs baseline-separated")
        
        print("\n" + "="*70)
        
        # Save to file if requested
        if save_path:
            with open(save_path, 'w') as f:
                f.write("QUALITY CONTROL REPORT\n")
                f.write("="*70 + "\n\n")
                
                if 'peak_shape' in self.qc_results:
                    f.write("Peak Shape Analysis:\n")
                    f.write(self.qc_results['peak_shape'].to_string(index=False))
                    f.write("\n\n")
                
                if 'resolution' in self.qc_results:
                    f.write("Peak Resolution Analysis:\n")
                    f.write(self.qc_results['resolution'].to_string(index=False))
                    f.write("\n")
            
            print(f"✅ QC report saved to: {save_path}")


# ============================================================================
# Complete Analysis Pipeline
# ============================================================================

class LCMSAnalysisPipeline:
    """
    Complete pipeline for LC-MS data analysis.
    """
    
    def __init__(self, output_dir='results'):
        """
        Initialize analysis pipeline.
        
        Args:
            output_dir: Directory for output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        print(f"\n{'='*70}")
        print(f"   LC-MS DATA PROCESSING PIPELINE")
        print(f"{'='*70}\n")
        print(f"📁 Output directory: {self.output_dir}")
    
    def process_sample(self, data_file, sample_name='Sample'):
        """
        Process a single LC-MS sample.
        
        Args:
            data_file: Path to chromatogram CSV file
            sample_name: Name for this sample
        """
        print(f"\n{'='*70}")
        print(f"Processing: {sample_name}")
        print(f"{'='*70}")
        
        results = {'sample_name': sample_name}
        
        # Parse chromatogram
        print("\n📊 Step 1: Loading Chromatogram")
        parser = ChromatogramParser(data_file)
        parser.plot_chromatogram(
            save_path=self.output_dir / f"{sample_name}_chromatogram.png"
        )
        results['chromatogram'] = parser.get_chromatogram()
        
        # Detect peaks
        print("\n📊 Step 2: Peak Detection")
        detector = PeakDetector(parser.get_chromatogram())
        detector.detect_peaks(prominence=0.1, width=5)
        detector.integrate_peaks()
        detector.calculate_snr()
        detector.plot_detected_peaks(
            save_path=self.output_dir / f"{sample_name}_peaks.png"
        )
        results['peaks'] = detector.peaks
        
        # Quality control
        print("\n📊 Step 3: Quality Control")
        qc = QualityControl()
        qc.check_peak_shape(detector.peaks, parser.get_chromatogram(), detector.baseline)
        qc.check_resolution(detector.peaks)
        qc.generate_qc_report(
            save_path=self.output_dir / f"{sample_name}_qc_report.txt"
        )
        results['qc'] = qc.qc_results
        
        # Save peak data
        detector.peaks.to_csv(
            self.output_dir / f"{sample_name}_peaks.csv", 
            index=False
        )
        
        print(f"\n{'='*70}")
        print(f"✅ {sample_name} processing complete!")
        print(f"{'='*70}\n")
        
        return results


# ============================================================================
# Example Usage
# ============================================================================

if __name__ == "__main__":
    print("""
    ╔═══════════════════════════════════════════════════════════════════╗
    ║     LC-MS DATA PROCESSING TOOLKIT                                 ║
    ╚═══════════════════════════════════════════════════════════════════╝
    
    This toolkit provides comprehensive LC-MS data analysis:
    
    1. Chromatogram parsing and visualization
    2. Automated peak detection and integration
    3. Calibration curve generation
    4. Sample quantification
    5. Quality control checks
    
    Example usage:
    
    # Complete pipeline
    pipeline = LCMSAnalysisPipeline(output_dir='results')
    results = pipeline.process_sample('data/sample1.csv', 'Sample_1')
    
    # Calibration and quantification
    cal = CalibrationCurve()
    cal.fit_calibration([1, 5, 10, 50, 100], [1000, 5000, 10000, 50000, 100000])
    cal.plot_calibration_curve()
    unknown_conc = cal.quantify_samples([7500, 25000])
    """)