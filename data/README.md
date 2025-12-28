# Data Directory

Place your LC-MS chromatogram CSV files here.

## Generating Example Data

Run the data generator script to create example chromatograms:

```bash
python generate_sample_data.py
```

This will create:
- `sample1.csv` - Good quality chromatogram with 4 peaks
- `sample2.csv` - Higher concentration sample
- `sample3.csv` - Noisy data example
- `standard_1ngml.csv` through `standard_100ngml.csv` - Calibration standards

## File Format

CSV files should have two columns:

```csv
time,intensity
0.00,1523.4
0.01,1589.2
0.02,1645.8
...
```

- **time**: Retention time in minutes (float)
- **intensity**: Signal intensity in arbitrary units (float)

## Exporting from LC-MS Software

### Thermo Xcalibur
1. Open chromatogram
2. File → Export → Select "Chromatogram"
3. Choose CSV format
4. Save with time and intensity columns

### Agilent MassHunter
1. Open chromatogram view
2. File → Export → Data
3. Select "Time vs. Intensity"
4. Export as CSV

### Waters Empower
1. Select chromatogram
2. File → Export
3. Choose "Comma Delimited" format

## Note

Large chromatogram files are not tracked by git. See `.gitignore`.