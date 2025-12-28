"""
Basic LC-MS Analysis Example

Demonstrates basic usage of the LC-MS Data Processor toolkit.
"""

import sys
sys.path.append('..')

from lcms_data_processor import LCMSAnalysisPipeline

print("""
╔═══════════════════════════════════════════════════════════════════╗
║     BASIC LC-MS ANALYSIS                                          ║
╚═══════════════════════════════════════════════════════════════════╝
""")

# Create pipeline
pipeline = LCMSAnalysisPipeline(output_dir='results')

# Process a sample
print("\nProcessing Sample 1...")
results = pipeline.process_sample(
    data_file='data/sample1.csv',
    sample_name='Sample_1'
)

# Display results summary
print("\n" + "="*70)
print("ANALYSIS SUMMARY")
print("="*70)

print(f"\n📊 Sample: {results['sample_name']}")
print(f"🔍 Peaks Detected: {len(results['peaks'])}")

if len(results['peaks']) > 0:
    print("\nPeak Details:")
    peak_summary = results['peaks'][['peak_id', 'retention_time', 'peak_area', 'snr']]
    print(peak_summary.to_string(index=False))
    
    # Highlight main peak
    main_peak = results['peaks'].loc[results['peaks']['peak_area'].idxmax()]
    print(f"\n⭐ Main Peak (highest area):")
    print(f"   RT: {main_peak['retention_time']:.2f} min")
    print(f"   Area: {main_peak['peak_area']:.2e}")
    print(f"   S/N: {main_peak['snr']:.1f}")

print("\n" + "="*70)
print("✅ Analysis Complete!")
print("📁 Check 'results/' folder for:")
print("   - Sample_1_chromatogram.png")
print("   - Sample_1_peaks.png")
print("   - Sample_1_peaks.csv")
print("   - Sample_1_qc_report.txt")
print("="*70 + "\n")