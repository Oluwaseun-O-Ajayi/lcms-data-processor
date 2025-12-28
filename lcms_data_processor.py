"""
LC-MS Data Processing Pipeline
===============================

A comprehensive Python toolkit for processing liquid chromatography-mass 
spectrometry (LC-MS) data. Designed for pharmaceutical, bioanalytical, and 
clinical chemistry laboratories.

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
import seaborn as sns
from scipy import signal, integrate
from scipy.optimize import curve_fit
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")


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
        
        plt.show()


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
    
    def detect_peaks(self, prominence=0.1, width=5, height=None):
        """
        Detect peaks using scipy's find_peaks.
        
        Args:
            prominence: Minimum peak prominence (relative to max intensity)
            width: Minimum peak width in data points
            height: Minimum peak height (absolute)
        """
        # Calculate baseline
        self.baseline = self._estimate_baseline()
        
        # Subtract baseline
        corrected_intensity = self.data['intensity'] - self.baseline
        
        # Set height threshold if not provided
        if height is None:
            height = corrected_intensity.max() * 0.05  # 5% of max
        
        # Detect peaks
        peak_indices, properties = signal.find_peaks(
            corrected_intensity,
            prominence=prominence * corrected_intensity.max(),
            width=width,
            height=height
        )
        
        # Extract peak information
        peaks_data = []
        for i, idx in enumerate(peak_indices):
            peak_info = {
                'peak_id': i + 1,
                'retention_time': self.data['time'].iloc[idx],
                'peak_height': corrected_intensity.iloc[idx],
                'prominence': properties['prominences'][i],
                'width': properties['widths'][i],
                'left_base': int(properties['left_bases'][i]),
                'right_base': int(properties['right_bases'][i])
            }
            peaks_data.append(peak_info)
        
        self.peaks = pd.DataFrame(peaks_data)
        
        return self.peaks
    
    def _estimate_baseline(self):
        """Estimate baseline using rolling minimum."""
        window = int(len(self.data) * 0.1)  # 10% of data points
        baseline = self.data['intensity'].rolling(
            window=window, center=True, min_periods=1
        ).min()
        return baseline
    
    def integrate_peaks(self):
        """Integrate detected peaks to get peak areas."""
        if self.peaks is None:
            raise ValueError("No peaks detected. Run detect_peaks() first.")
        
        areas = []
        
        for _, peak in self.peaks.iterrows():
            # Get peak region
            left = peak['left_base']
            right = peak['right_base']
            
            peak_time = self.data['time'].iloc[left:right+1]
            peak_intensity = (self.data['intensity'].iloc[left:right+1] - 
                            self.baseline.iloc[left:right+1])
            
            # Integrate using trapezoidal rule
            area = integrate.trapezoid(peak_intensity, peak_time)
            areas.append(area)
        
        self.peaks['peak_area'] = areas
        
        return self.peaks
    
    def calculate_snr(self):
        """Calculate signal-to-noise ratio for each peak."""
        if self.peaks is None:
            raise ValueError("No peaks detected. Run detect_peaks() first.")
        
        # Estimate noise from baseline region
        noise_region = self.data['intensity'][:100]  # First 100 points
        noise_std = noise_region.std()
        
        self.peaks['snr'] = self.peaks['peak_height'] / noise_std
        
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
        
        ax2.set_xlabel('Peak ID', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Peak Area', fontsize=11, fontweight='bold', color='skyblue')
        ax2_twin.set_ylabel('S/N Ratio', fontsize=11, fontweight='bold', color='coral')
        ax2.set_title('Peak Metrics', fontsize=13, fontweight='bold')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(self.peaks['peak_id'])
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Combine legends
        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2_twin.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✅ Peak detection plot saved to: {save_path}")
        
        plt.show()
        
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
        
        # Fit linear model: y = mx + b
        self.curve_params = np.polyfit(concentrations, peak_areas, 1)
        
        # Calculate R²
        predicted_areas = np.polyval(self.curve_params, concentrations)
        ss_res = np.sum((peak_areas - predicted_areas) ** 2)
        ss_tot = np.sum((peak_areas - np.mean(peak_areas)) ** 2)
        self.r_squared = 1 - (ss_res / ss_tot)
        
        return self.curve_params, self.r_squared
    
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
        
        plt.show()
        
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
    
    def check_peak_shape(self, peak_data, chromatogram_data):
        """
        Assess peak shape quality (symmetry).
        
        Args:
            peak_data: DataFrame of detected peaks
            chromatogram_data: Full chromatogram DataFrame
        """
        shape_metrics = []
        
        for _, peak in peak_data.iterrows():
            # Get peak region
            left = peak['left_base']
            right = peak['right_base']
            peak_idx = int((left + right) / 2)
            
            # Calculate asymmetry factor
            # AF = b/a where b is distance from peak max to trailing edge
            # and a is distance from leading edge to peak max
            a = peak_idx - left
            b = right - peak_idx
            asymmetry = b / a if a > 0 else np.nan
            
            shape_metrics.append({
                'peak_id': peak['peak_id'],
                'asymmetry_factor': asymmetry,
                'quality': 'Good' if 0.9 <= asymmetry <= 1.1 else 
                          'Acceptable' if 0.8 <= asymmetry <= 1.2 else 'Poor'
            })
        
        self.qc_results['peak_shape'] = pd.DataFrame(shape_metrics)
        return self.qc_results['peak_shape']
    
    def check_resolution(self, peak_data):
        """
        Calculate resolution between adjacent peaks.
        
        Args:
            peak_data: DataFrame of detected peaks
        """
        if len(peak_data) < 2:
            print("⚠️ Need at least 2 peaks to calculate resolution")
            return None
        
        resolutions = []
        
        for i in range(len(peak_data) - 1):
            peak1 = peak_data.iloc[i]
            peak2 = peak_data.iloc[i + 1]
            
            # Resolution = 2 * (RT2 - RT1) / (W1 + W2)
            rt_diff = peak2['retention_time'] - peak1['retention_time']
            avg_width = (peak1['width'] + peak2['width']) / 2
            
            resolution = 2 * rt_diff / (avg_width * 0.1) if avg_width > 0 else np.nan
            
            resolutions.append({
                'peak_pair': f"{peak1['peak_id']}-{peak2['peak_id']}",
                'resolution': resolution,
                'quality': 'Baseline' if resolution >= 1.5 else
                          'Partial' if resolution >= 1.0 else 'Poor'
            })
        
        self.qc_results['resolution'] = pd.DataFrame(resolutions)
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
        qc.check_peak_shape(detector.peaks, parser.get_chromatogram())
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