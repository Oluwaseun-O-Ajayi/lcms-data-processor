# LC-MS Data Processing Pipeline

Open-source Python toolkit for reproducible LC-MS chromatogram processing, peak integration, calibration, quantification, and quality-control analysis

The framework provides chromatogram visualization, baseline correction, peak detection, numerical integration, calibration-curve modeling, sample quantification, quality-control reporting, and publication-quality figure generation for analytical chemistry, bioanalysis, and laboratory research workflows

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

LC-MS Data Processor is a reproducible computational workflow for chromatogram analysis using exported time–intensity data

The software integrates:

- **Chromatogram visualization**
- **Baseline estimation**
- **Automated peak detection**
- **Peak integration**
- **Signal-to-noise calculation**
- **Calibration-curve modeling**
- **Sample quantification**
- **Quality-control assessment**
- **Automated reporting**

The framework is intended for educational, academic, industrial, and analytical chemistry environments where transparent and reproducible chromatogram-processing workflows are required

## Scope and Intended Use

This repository provides computational workflows for LC-MS chromatogram processing and quantitative data analysis

Implemented capabilities include:

- **Peak detection**
- **Peak integration**
- **Calibration modeling**
- **Quantification**
- **Quality-control assessment**
- **Reporting**

The software operates on exported chromatogram data and is intended for research, education, workflow development, and reproducibility-focused applications

Users are responsible for validating analytical methods, calibration strategies, integration settings, and quality-control criteria before applying results to regulated or decision-critical workflows

## Scientific Software Contribution

The primary contribution of this repository is a reusable software workflow for chromatogram-based LC-MS analysis

Implemented capabilities include:

- **Chromatogram parsing**
- **Baseline correction**
- **Automated peak detection**
- **Trapezoidal peak integration**
- **Signal-to-noise analysis**
- **Calibration-curve fitting**
- **Concentration prediction**
- **Peak-shape assessment**
- **Peak-resolution assessment**
- **Automated report generation**

The framework is designed to improve accessibility and reproducibility of LC-MS data analysis workflows

## Features

-  **Automated peak detection** - Uses scipy's signal processing algorithms
-  **Baseline correction** - Rolling minimum baseline estimation
-  **Peak integration** - Trapezoidal integration for accurate peak areas
-  **Calibration curves** - Linear regression with goodness-of-fit metrics
-  **Quality metrics** - S/N ratio, peak asymmetry, resolution
-  **Batch processing** - Analyze multiple samples in sequence
-  **Report generation** - Automated QC reports and summary tables
-  **Publication-quality plots** - High-resolution figures for papers

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/Oluwaseun-O-Ajayi/lcms-data-processor.git
cd lcms-data-processor

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```python
from lcms_data_processor import LCMSAnalysisPipeline

# Process a single sample
pipeline = LCMSAnalysisPipeline(output_dir='results')
results = pipeline.process_sample('data/sample1.csv', sample_name='Sample_1')
```

This will:
- Load and plot the chromatogram
- Detect and integrate peaks
- Calculate S/N ratios
- Perform quality control checks
- Generate all plots and reports

## Requirements

- Python 3.8+
- numpy >= 1.23.0
- pandas >= 1.5.0
- matplotlib >= 3.6.0
- seaborn >= 0.12.0
- scipy >= 1.9.0

## Data Format

Input CSV files should have two columns:

```csv
time,intensity
0.00,1523.4
0.01,1589.2
0.02,1645.8
...
```

- **time**: Retention time in minutes
- **intensity**: Signal intensity in arbitrary units (AU)

## Example Data

The repository includes synthetic chromatogram datasets generated from Gaussian peaks, baseline signal, and controlled noise

These datasets are provided for:

- **Software testing**
- **Workflow demonstrations**
- **Educational purposes**
- **Reproducibility examples**

Example chromatograms should not be interpreted as validated experimental LC-MS measurements

## Available Modules

### 1. ChromatogramParser

Load and visualize chromatogram data.

```python
from lcms_data_processor import ChromatogramParser

# Load chromatogram
parser = ChromatogramParser('data/sample.csv')

# Get summary
metadata = parser.metadata
print(f"Time range: {metadata['time_range']}")
print(f"Max intensity: {metadata['max_intensity']}")

# Plot chromatogram
parser.plot_chromatogram(save_path='chromatogram.png')

# Get data for custom analysis
df = parser.get_chromatogram()
```

### 2. PeakDetector

Detect, integrate, and characterize peaks.

```python
from lcms_data_processor import ChromatogramParser, PeakDetector

parser = ChromatogramParser('data/sample.csv')
detector = PeakDetector(parser.get_chromatogram())

# Detect peaks
peaks = detector.detect_peaks(
    prominence=0.1,  # 10% of max intensity
    width=5,         # Minimum width in data points
    height=None      # Auto-calculate
)

# Integrate peaks
peaks_with_areas = detector.integrate_peaks()

# Calculate S/N ratios
peaks_with_snr = detector.calculate_snr()

# Visualize
detector.plot_detected_peaks(save_path='peaks.png')

# Export peak data
peaks.to_csv('peak_data.csv', index=False)
```

**Output metrics:**
- Peak ID
- Retention time
- Peak height
- Peak area
- Signal-to-noise ratio (S/N)
- Peak width

### 3. CalibrationCurve

Generate calibration curves and quantify unknowns.

```python
from lcms_data_processor import CalibrationCurve

cal = CalibrationCurve()

# Fit calibration curve with standards
concentrations = [1, 5, 10, 50, 100]  # ng/mL
peak_areas = [1000, 5000, 10000, 50000, 100000]

slope, r_squared = cal.fit_calibration(concentrations, peak_areas)

# Visualize calibration
cal.plot_calibration_curve(save_path='calibration.png')

# Quantify unknown samples
unknown_areas = [7500, 25000, 75000]
unknown_concentrations = cal.quantify_samples(unknown_areas)

print(f"Unknown concentrations: {unknown_concentrations}")
```

**Output:**
- Slope and intercept
- R² value
- Concentration predictions

### 4. QualityControl

Assess data quality with industry-standard metrics.

```python
from lcms_data_processor import QualityControl

qc = QualityControl()

# Check peak shape (asymmetry factor)
shape_qc = qc.check_peak_shape(peaks_df, chromatogram_df)

# Check peak resolution
resolution_qc = qc.check_resolution(peaks_df)

# Generate QC report
qc.generate_qc_report(save_path='qc_report.txt')
```

**Metrics:**
- **Peak asymmetry** - Measures peak tailing (ideal: 0.9-1.1)
- **Resolution** - Separation between adjacent peaks (ideal: >1.5)
- **S/N ratio** - Signal quality (ideal: >10)

### 5. Complete Pipeline

Run all analyses at once.

```python
from lcms_data_processor import LCMSAnalysisPipeline

pipeline = LCMSAnalysisPipeline(output_dir='results')

# Process single sample
results = pipeline.process_sample(
    data_file='data/sample1.csv',
    sample_name='Control_1'
)

# Access results
print(f"Detected {len(results['peaks'])} peaks")
print(results['peaks'])
```

##  Project Structure

```
lcms-data-processor/
├── lcms_data_processor.py     # Main toolkit
├── README.md                   # This file
├── requirements.txt            # Dependencies
├── LICENSE                     # MIT License
├── .gitignore                 # Git ignore rules
├── examples/                   # Example scripts
│   ├── basic_analysis.py
│   ├── calibration_workflow.py
│   └── batch_processing.py
├── data/                       # Example data
│   ├── sample1.csv
│   ├── sample2.csv
│   └── README.md
└── results/                    # Output directory
    └── README.md
```

## Example Workflows

### Workflow 1: Single Sample Analysis

```python
from lcms_data_processor import LCMSAnalysisPipeline

# Quick analysis of one sample
pipeline = LCMSAnalysisPipeline()
results = pipeline.process_sample('data/sample1.csv', 'Test_Sample')

# Results include:
# - Chromatogram plot
# - Peak detection plot
# - QC report
# - Peak data CSV
```

### Workflow 2: Calibration and Quantification

```python
from lcms_data_processor import CalibrationCurve, PeakDetector, ChromatogramParser

# 1. Build calibration curve from standards
cal = CalibrationCurve()

standards_conc = [1, 5, 10, 25, 50, 100]  # ng/mL
standards_areas = []

for std_file in standard_files:
    parser = ChromatogramParser(std_file)
    detector = PeakDetector(parser.get_chromatogram())
    peaks = detector.detect_peaks()
    detector.integrate_peaks()
    
    # Get area of analyte peak (e.g., peak at RT ~5.2 min)
    analyte_peak = peaks[peaks['retention_time'].between(5.0, 5.4)]
    standards_areas.append(analyte_peak['peak_area'].values[0])

cal.fit_calibration(standards_conc, standards_areas)
cal.plot_calibration_curve(save_path='calibration.png')

# 2. Quantify unknowns
for unknown_file in unknown_files:
    parser = ChromatogramParser(unknown_file)
    detector = PeakDetector(parser.get_chromatogram())
    peaks = detector.detect_peaks()
    detector.integrate_peaks()
    
    analyte_peak = peaks[peaks['retention_time'].between(5.0, 5.4)]
    peak_area = analyte_peak['peak_area'].values[0]
    
    concentration = cal.quantify_samples([peak_area])[0]
    print(f"{unknown_file}: {concentration:.2f} ng/mL")
```

### Workflow 3: Batch Processing

```python
from lcms_data_processor import LCMSAnalysisPipeline
import glob

pipeline = LCMSAnalysisPipeline(output_dir='batch_results')

# Process all CSV files in data folder
data_files = glob.glob('data/*.csv')

batch_results = []
for data_file in data_files:
    sample_name = Path(data_file).stem
    results = pipeline.process_sample(data_file, sample_name)
    
    # Collect summary data
    batch_results.append({
        'sample': sample_name,
        'num_peaks': len(results['peaks']),
        'total_peak_area': results['peaks']['peak_area'].sum(),
        'avg_snr': results['peaks']['snr'].mean()
    })

# Create summary report
import pandas as pd
summary = pd.DataFrame(batch_results)
summary.to_csv('batch_results/summary.csv', index=False)
print(summary)
```

## Research Applications
The examples included in this repository demonstrate workflow capabilities and are not intended to substitute for laboratory-specific method validation requirements

### Pharmaceutical Analysis
- **Drug quantification** - Measure API concentrations in formulations
- **Impurity profiling** - Detect and quantify degradation products
- **Stability studies** - Monitor drug degradation over time
- **Quality control** - Verify batch-to-batch consistency

### Bioanalytical Chemistry
- **Pharmacokinetics** - Quantify drug levels in plasma/serum
- **Metabolite identification** - Detect drug metabolites
- **Biomarker analysis** - Measure disease biomarkers
- **Therapeutic drug monitoring** - Ensure proper dosing

### Clinical Chemistry
- **Small molecule quantification** - Hormones, vitamins, toxins
- **Screening assays** - Newborn screening, drug testing
- **Method validation** - Linearity, accuracy, precision

### Environmental Analysis
- **Contaminant detection** - Pesticides, pollutants
- **Water quality testing** - PFAS, pharmaceuticals in water
- **Food safety** - Residue analysis

## Customization

### Adjust Peak Detection Sensitivity

```python
# More sensitive (detect smaller peaks)
detector.detect_peaks(prominence=0.05, width=3, height=100)

# Less sensitive (only major peaks)
detector.detect_peaks(prominence=0.2, width=10, height=1000)
```

### Custom Baseline Correction

```python
# Modify baseline window size
# In PeakDetector._estimate_baseline():
window = int(len(self.data) * 0.05)  # 5% instead of 10%
```

### High-Resolution Plots

```python
# For publication
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 600

detector.plot_detected_peaks(save_path='high_res_peaks.png')
```

## Quality Metrics Guide

### Signal-to-Noise Ratio (S/N)
- **S/N > 10**: Excellent signal
- **S/N 3-10**: Acceptable
- **S/N < 3**: Poor, may need method optimization

### Peak Asymmetry Factor
- **0.9-1.1**: Symmetric peak (ideal)
- **0.8-1.2**: Acceptable
- **<0.8 or >1.2**: Tailing or fronting (investigate column/method)

### Resolution
- **Rs ≥ 1.5**: Baseline separation (ideal)
- **Rs 1.0-1.5**: Partial separation
- **Rs < 1.0**: Poor separation (co-elution risk)

## Troubleshooting

### No Peaks Detected
- Lower `prominence` parameter
- Check if data loaded correctly
- Verify intensity units are consistent

### Too Many False Peaks
- Increase `prominence` and `width` parameters
- Improve baseline correction
- Filter data before analysis

### Poor Calibration (R² < 0.95)
- Check for outliers in standards
- Verify concentration range is appropriate
- Ensure peak integration is consistent

## Integration with Other Tools

Works well with:
- **Thermo Xcalibur** - Export chromatograms as CSV
- **Agilent MassHunter** - Export time vs. intensity data
- **Waters Empower** - Export results to CSV format
- **Python ecosystem** - Combine with pandas, scikit-learn

## Validation and Interpretation

Model outputs should be interpreted alongside analytical method validation and experimental quality considerations

Recommended practices include:

- **Visual inspection of chromatograms**
- **Review of integration boundaries**
- **Verification of calibration performance**
- **Assessment of signal-to-noise ratios**
- **Evaluation of peak-shape metrics**
- **Confirmation of peak resolution**

Software outputs should support, not replace, analytical judgment

## Contributing

Contributions welcome! Areas for enhancement:
- Advanced baseline algorithms (e.g., asymmetric least squares)
- Non-linear calibration models (quadratic, power)
- Peak deconvolution for overlapping peaks
- Integration with vendor software APIs
- Support for additional file formats (mzML, mzXML)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

##  Acknowledgments

- Developed for pharmaceutical and bioanalytical chemistry research
- Inspired by industry-standard HPLC/LC-MS software
- Thanks to scipy and the scientific Python community

## Contact

**Oluwaseun O. Ajayi**  
PhD Researcher, Chemistry  
University of Georgia

- **GitHub**: [@Oluwaseun-O-Ajayi](https://github.com/Oluwaseun-O-Ajayi)
- **Academic Email**: oluwaseun.ajayi@uga.edu
- **Personal Email**: seunolanikeajayi@gmail.com

## Publication

Associated manuscript:

Ajayi OO.

LC-MS Data Processor: An Open-Source Python Workflow for Chromatogram Processing, Peak Integration, Calibration, and Quality-Control Reporting

Manuscript in preparation

## Citation

If you use this toolkit in your research:

```bibtex
@software{lcms_data_processor,
  author = {Oluwaseun O. Ajayi},
  title = {LC-MS Data Processing Pipeline},
  year = {2024},
  url = {https://github.com/Oluwaseun-O-Ajayi/lcms-data-processor}
}
```

---

**Advancing reproducible LC-MS data analysis through open scientific software**
