# Results Directory

Analysis outputs will be saved here.

## Generated Files

When you run the analysis scripts, the following files will be created:

### Chromatogram Plots
- `{sample_name}_chromatogram.png` - Raw chromatogram with metadata
- `{sample_name}_peaks.png` - Detected peaks with integration regions
- `calibration_curve.png` - Calibration curve with standards

### Data Files
- `{sample_name}_peaks.csv` - Peak data table (RT, area, S/N, etc.)
- `{sample_name}_qc_report.txt` - Quality control metrics

### Batch Analysis
- `batch_summary.csv` - Summary table for multiple samples

## Running Examples

Generate example data first:
```bash
python generate_sample_data.py
```

Then run analyses:
```bash
# Basic analysis
python -m examples.basic_analysis

# Calibration workflow
python -m examples.calibration_workflow

# Batch processing
python -m examples.batch_processing
```

## Note

Generated output files (PNG, CSV, TXT) are not tracked by git to keep the repository size small. See `.gitignore`.