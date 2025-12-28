"""
Calibration and Quantification Workflow

Demonstrates how to build a calibration curve and quantify unknowns.
"""

import sys
sys.path.append('..')

from lcms_data_processor import (
    ChromatogramParser,
    PeakDetector,
    CalibrationCurve
)
import glob

print("""
╔═══════════════════════════════════════════════════════════════════╗
║     CALIBRATION & QUANTIFICATION WORKFLOW                         ║
╚═══════════════════════════════════════════════════════════════════╝
""")

# Step 1: Process calibration standards
print("\n" + "="*70)
print("STEP 1: Processing Calibration Standards")
print("="*70)

standard_files = sorted(glob.glob('data/standard_*.csv'))

if len(standard_files) == 0:
    print("\n⚠️  No standard files found!")
    print("Run 'python generate_sample_data.py' first to create example data.")
    sys.exit(1)

standards_data = []

for std_file in standard_files:
    # Extract concentration from filename (e.g., "standard_50ngml.csv" -> 50)
    conc_str = std_file.split('_')[1].replace('ngml.csv', '')
    concentration = float(conc_str)
    
    # Parse and analyze
    parser = ChromatogramParser(std_file)
    detector = PeakDetector(parser.get_chromatogram())
    detector.detect_peaks(prominence=0.1, width=5)
    detector.integrate_peaks()
    
    # Get analyte peak (should be around RT 5.2 min)
    analyte_peaks = detector.peaks[
        detector.peaks['retention_time'].between(5.0, 5.5)
    ]
    
    if len(analyte_peaks) > 0:
        peak_area = analyte_peaks['peak_area'].values[0]
        standards_data.append({
            'concentration': concentration,
            'peak_area': peak_area
        })
        print(f"✅ {std_file}: {concentration} ng/mL -> Area = {peak_area:.2e}")
    else:
        print(f"⚠️  {std_file}: No analyte peak found")

# Step 2: Build calibration curve
print("\n" + "="*70)
print("STEP 2: Building Calibration Curve")
print("="*70)

if len(standards_data) < 3:
    print("\n⚠️  Need at least 3 standards for calibration")
    sys.exit(1)

cal = CalibrationCurve()

concentrations = [s['concentration'] for s in standards_data]
peak_areas = [s['peak_area'] for s in standards_data]

slope, r_squared = cal.fit_calibration(concentrations, peak_areas)

print(f"\n📊 Calibration Results:")
print(f"   Slope: {slope:.4e}")
print(f"   R²: {r_squared:.4f}")

if r_squared >= 0.99:
    print("   ✅ Excellent calibration!")
elif r_squared >= 0.95:
    print("   ✅ Good calibration")
else:
    print("   ⚠️  Warning: R² < 0.95")

cal.plot_calibration_curve(save_path='results/calibration_curve.png')

# Step 3: Quantify unknown samples
print("\n" + "="*70)
print("STEP 3: Quantifying Unknown Samples")
print("="*70)

unknown_files = ['data/sample1.csv', 'data/sample2.csv']

for unknown_file in unknown_files:
    try:
        # Parse and analyze
        parser = ChromatogramParser(unknown_file)
        detector = PeakDetector(parser.get_chromatogram())
        detector.detect_peaks(prominence=0.1, width=5)
        detector.integrate_peaks()
        
        # Get analyte peak
        analyte_peaks = detector.peaks[
            detector.peaks['retention_time'].between(5.0, 5.5)
        ]
        
        if len(analyte_peaks) > 0:
            peak_area = analyte_peaks['peak_area'].values[0]
            concentration = cal.quantify_samples([peak_area])[0]
            
            print(f"\n📊 {unknown_file}:")
            print(f"   Peak Area: {peak_area:.2e}")
            print(f"   Concentration: {concentration:.2f} ng/mL")
        else:
            print(f"\n⚠️  {unknown_file}: No analyte peak found")
            
    except FileNotFoundError:
        print(f"\n⚠️  {unknown_file}: File not found")
        print("   Run 'python generate_sample_data.py' to create sample data")

print("\n" + "="*70)
print("✅ Calibration Workflow Complete!")
print("📁 Check 'results/calibration_curve.png'")
print("="*70 + "\n")