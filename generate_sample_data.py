"""
Generate Example LC-MS Chromatogram Data

Creates synthetic chromatogram data for testing the LC-MS toolkit.
"""

import numpy as np
import pandas as pd
from pathlib import Path

def gaussian_peak(x, amplitude, center, width):
    """Generate a Gaussian peak."""
    return amplitude * np.exp(-((x - center) ** 2) / (2 * width ** 2))

def generate_chromatogram(
    time_range=(0, 15),
    num_points=1500,
    peaks=None,
    baseline=1000,
    noise_level=50
):
    """
    Generate synthetic chromatogram data.
    
    Args:
        time_range: (start, end) time in minutes
        num_points: Number of data points
        peaks: List of (retention_time, height, width) tuples
        baseline: Baseline intensity level
        noise_level: Standard deviation of noise
    
    Returns:
        DataFrame with 'time' and 'intensity' columns
    """
    if peaks is None:
        # Default peaks
        peaks = [
            (3.5, 5000, 0.15),   # Small peak
            (5.2, 15000, 0.2),   # Main peak (analyte)
            (7.8, 8000, 0.25),   # Medium peak
            (10.5, 3000, 0.18),  # Small peak
        ]
    
    # Generate time points
    time = np.linspace(time_range[0], time_range[1], num_points)
    
    # Start with baseline
    intensity = np.ones(num_points) * baseline
    
    # Add peaks
    for rt, height, width in peaks:
        intensity += gaussian_peak(time, height, rt, width)
    
    # Add noise
    noise = np.random.normal(0, noise_level, num_points)
    intensity += noise
    
    # Ensure no negative values
    intensity = np.maximum(intensity, 0)
    
    # Create DataFrame
    df = pd.DataFrame({
        'time': time,
        'intensity': intensity
    })
    
    return df

def main():
    """Generate example chromatogram files."""
    
    # Create data directory
    data_dir = Path('data')
    data_dir.mkdir(exist_ok=True)
    
    print("Generating example LC-MS chromatogram data...")
    
    # Sample 1: Good quality chromatogram
    print("\n1. Generating sample1.csv (good quality)...")
    chrom1 = generate_chromatogram(
        peaks=[
            (3.5, 5000, 0.15),
            (5.2, 15000, 0.2),
            (7.8, 8000, 0.25),
            (10.5, 3000, 0.18),
        ],
        baseline=1000,
        noise_level=50
    )
    chrom1.to_csv(data_dir / 'sample1.csv', index=False)
    print(f"   ✅ Created: {data_dir / 'sample1.csv'}")
    print(f"   Peaks: 4, S/N: High")
    
    # Sample 2: Higher concentration (for calibration)
    print("\n2. Generating sample2.csv (higher concentration)...")
    chrom2 = generate_chromatogram(
        peaks=[
            (3.5, 7000, 0.15),
            (5.2, 25000, 0.2),   # Higher main peak
            (7.8, 12000, 0.25),
            (10.5, 4500, 0.18),
        ],
        baseline=1000,
        noise_level=50
    )
    chrom2.to_csv(data_dir / 'sample2.csv', index=False)
    print(f"   ✅ Created: {data_dir / 'sample2.csv'}")
    print(f"   Peaks: 4, S/N: Very High")
    
    # Sample 3: Lower quality (more noise)
    print("\n3. Generating sample3.csv (noisy data)...")
    chrom3 = generate_chromatogram(
        peaks=[
            (3.5, 4000, 0.15),
            (5.2, 10000, 0.2),
            (7.8, 6000, 0.25),
        ],
        baseline=1000,
        noise_level=150  # Higher noise
    )
    chrom3.to_csv(data_dir / 'sample3.csv', index=False)
    print(f"   ✅ Created: {data_dir / 'sample3.csv'}")
    print(f"   Peaks: 3, S/N: Medium (noisy)")
    
    # Calibration standards
    print("\n4. Generating calibration standards...")
    standard_concentrations = [1, 5, 10, 50, 100]  # ng/mL
    
    for i, conc in enumerate(standard_concentrations, 1):
        # Scale peak height proportionally to concentration
        main_peak_height = conc * 150
        
        chrom = generate_chromatogram(
            peaks=[
                (5.2, main_peak_height, 0.2),  # Only analyte peak
            ],
            baseline=1000,
            noise_level=30
        )
        
        filename = f'standard_{conc}ngml.csv'
        chrom.to_csv(data_dir / filename, index=False)
        print(f"   ✅ Created: {data_dir / filename} ({conc} ng/mL)")
    
    print("\n" + "="*60)
    print("✅ Sample data generation complete!")
    print(f"📁 Files saved in: {data_dir.absolute()}")
    print("="*60)
    print("\nYou can now run:")
    print("  python -m examples.basic_analysis")
    print("  python -m examples.calibration_workflow")
    print("  python -m examples.batch_processing")

if __name__ == "__main__":
    main()