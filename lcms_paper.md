---
title: 'LC-MS Data Processor: An Open-Source Python Workflow for Chromatogram Processing, Peak Integration, Calibration, and Quality-Control Reporting'
tags:
  - Python
  - LC-MS
  - liquid chromatography
  - mass spectrometry
  - chromatogram processing
  - peak detection
  - analytical chemistry
  - bioanalytical chemistry
  - reproducible research
authors:
  - name: Oluwaseun O. Ajayi
    orcid: 0000-0003-0040-7217
    affiliation: 1
affiliations:
 - name: Department of Chemistry, University of Georgia, Athens, Georgia, USA
   index: 1
date: 2026-07-28
bibliography: paper.bib
---

# Summary

Liquid chromatography-mass spectrometry (LC-MS) is widely used in analytical chemistry, pharmaceutical analysis, bioanalytical chemistry, clinical chemistry, environmental monitoring, and related fields. Routine LC-MS data processing commonly involves chromatogram visualization, peak detection, peak integration, calibration-curve construction, analyte quantification, quality-control review, and report preparation. In many educational and resource-limited research settings, these tasks may be distributed across vendor software, spreadsheet calculations, custom scripts, and manually assembled reports. This fragmented workflow can limit reproducibility, reduce transparency, and make it difficult for students or small research groups to inspect and adapt the analysis logic.

`LC-MS Data Processor` is an open-source Python toolkit that provides a reproducible workflow for LC-MS chromatogram analysis using simple time-intensity CSV files. The software includes modules for chromatogram parsing, baseline estimation, automated peak detection, trapezoidal peak integration, signal-to-noise calculation, linear calibration-curve fitting, analyte quantification from peak area, peak-shape assessment, peak-resolution assessment, batch-style processing, and automated output generation. The package is designed to generate readable figures, CSV outputs, and quality-control reports from a script-driven workflow.

The current release is intended as an accessible framework for teaching, prototyping, method-development demonstrations, and reproducibility-focused LC-MS data analysis. Example chromatograms are generated synthetically using Gaussian peaks, baseline signal, and random noise to demonstrate expected input formats and processing behavior. These demonstration datasets are not presented as validated experimental LC-MS measurements. Instead, they support software testing, educational use, and transparent documentation of the computational workflow.

# Statement of need

LC-MS analysis is a central tool in quantitative analytical chemistry, but practical data-processing workflows can be difficult to reproduce outside vendor-specific software environments. Many LC-MS instruments are paired with proprietary platforms that are powerful but may not be available away from instrument workstations or to every student, collaborator, or resource-limited laboratory. Researchers may therefore rely on exported chromatogram tables, spreadsheets, screenshots, or separate scripts for downstream review and reporting. These approaches can support analysis, but they also increase the risk of inconsistent processing choices and manual transcription errors.

Open, inspectable workflows are especially valuable for students learning chromatogram analysis, researchers developing small-scale assays, and laboratories that need transparent computational steps before adopting larger informatics systems. A lightweight Python framework can help users understand how chromatogram data move from raw time-intensity traces to detected peaks, integrated areas, calibration curves, reported concentrations, and quality-control summaries.

`LC-MS Data Processor` addresses this need by collecting common LC-MS processing operations into a single, documented Python workflow. The toolkit reads CSV chromatogram files containing retention time and intensity columns; estimates a rolling-minimum baseline; applies signal-based peak detection; integrates peak areas using the trapezoidal rule; calculates signal-to-noise estimates; fits linear calibration curves; quantifies unknown samples from peak areas; computes peak-shape and resolution checks; and writes plots, tables, and textual quality-control reports. The contribution is not a new LC-MS theory or a replacement for validated vendor workflows, but an accessible and reproducible software implementation for education, prototyping, and transparent analytical workflow development.

# Software functionality

The package currently provides the following workflow components:

- Chromatogram parsing from CSV files with retention time and intensity columns.
- Metadata extraction, including filename, number of points, time range, and maximum intensity.
- Raw chromatogram visualization with annotated metadata.
- Rolling-minimum baseline estimation and baseline-corrected intensity traces.
- Automated peak detection using signal-processing routines.
- Peak integration using trapezoidal numerical integration.
- Signal-to-noise calculation based on baseline-region noise estimates.
- Peak visualization with retention-time markers, integration information, peak areas, and signal-to-noise summaries.
- Linear calibration-curve fitting from standard concentrations and peak areas.
- R-squared calculation for calibration goodness-of-fit.
- Quantification of unknown peak areas using the fitted calibration model.
- Peak-shape and adjacent-peak resolution checks.
- Single-sample pipeline processing with plot, CSV, and quality-control report outputs.
- Example scripts for basic analysis, calibration workflows, and batch-style processing.

# Demonstration data and example workflows

The repository includes a synthetic data-generation script that creates example chromatograms using Gaussian peaks, baseline intensity, and random noise. The script produces good-quality, higher-concentration, and noisier chromatogram examples, along with calibration-standard examples across a concentration series. These files are intended to demonstrate input formatting, processing behavior, and output generation.

Example workflows currently include single-sample analysis, calibration and quantification, and batch-style processing. In the single-sample workflow, a chromatogram is loaded, plotted, processed for peak detection, integrated, reviewed for signal-to-noise, assessed using quality-control checks, and exported as figures, peak tables, and a text report. In the calibration workflow, standard chromatograms are processed, analyte peak areas are extracted near a target retention-time window, a linear calibration curve is fitted, and unknown sample concentrations are calculated from integrated peak areas.

# Research and teaching applications

`LC-MS Data Processor` is most appropriate for educational, prototyping, and reproducibility-focused workflows in analytical chemistry. Potential use cases include teaching chromatogram processing concepts, demonstrating peak detection and integration, prototyping calibration workflows from exported CSV files, comparing basic quality-control metrics across samples, and creating transparent examples for students or early-stage researchers who need to understand how quantitative LC-MS outputs are produced.

The software may also provide a foundation for future extensions toward more advanced analytical workflows, including alternative baseline-correction methods, non-linear calibration models, overlapping-peak deconvolution, additional file-format support, and integration with open mass-spectrometry data standards. Such extensions should be validated carefully before use in regulated or decision-critical settings.

# Limitations

The current version is a lightweight open-source workflow and should not be interpreted as a validated replacement for regulated LC-MS processing systems or vendor-validated quantitative software. The included example chromatograms are synthetic demonstration datasets rather than experimental LC-MS measurements. Users are responsible for validating peak-detection settings, baseline-correction assumptions, integration boundaries, calibration-model suitability, signal-to-noise definitions, and quality-control thresholds for their own analytical methods. The package currently operates on exported CSV time-intensity data rather than native vendor files or open mass-spectrometry formats such as mzML.

# Acknowledgements

The author acknowledges the scientific Python ecosystem, including NumPy, SciPy, pandas, Matplotlib, and seaborn, which provide foundational tools for numerical computation, signal processing, data handling, and visualization. The development of this software was motivated by analytical chemistry and bioanalytical workflow needs, particularly the need for accessible and transparent computational tools for students, academic laboratories, and resource-limited research environments.

# References
