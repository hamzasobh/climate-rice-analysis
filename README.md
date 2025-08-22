# Climate Change Impact on Rice Production in Kafr Elsheikh

This repository contains Python code for analyzing climate data and its impact on rice production in the Kafr Elsheikh region of Egypt.

## Project Overview

This project analyzes:
- Historical and projected temperature patterns
- Precipitation variability
- Water vapor pressure calculations
- Regional climate averages for the Nile Delta
- Climate-agriculture relationships

## Repository Structure

```
climate-rice-analysis/
├── src/                    # Source code
│   └── climate_data_analysis.py
├── data/                   # Data files (add your .tif files here)
├── results/                # Output plots and analysis results
├── docs/                   # Documentation
├── requirements.txt        # Python dependencies
└── README.md              # This file
```

## Features

### Climate Data Analysis Functions

1. **Temperature Analysis**
   - `plot_temperature_map()` - Visualize temperature patterns
   - `plot_regional_temperature_analysis()` - Regional temperature analysis with averages
   - `extract_point_value()` - Extract values at specific coordinates

2. **Precipitation Analysis**
   - `plot_precipitation_map()` - Visualize precipitation patterns
   - Regional precipitation statistics

3. **Vapor Pressure Analysis**
   - `saturated_vapor_pressure()` - Calculate vapor pressure using Magnus-Tetens formula
   - `analyze_vapor_pressure()` - Analyze and visualize vapor pressure

4. **Regional Analysis**
   - `calculate_regional_average()` - Calculate averages for defined regions
   - Nile Delta specific analysis functions

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/climate-rice-analysis.git
cd climate-rice-analysis
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```python
from src.climate_data_analysis import *

# Analyze temperature data
plot_temperature_map("data/temperature_file.tif", "Temperature Analysis", "results/temp_map.png")

# Analyze precipitation
plot_precipitation_map("data/precipitation_file.tif", "Precipitation Analysis", "results/precip_map.png")

# Calculate regional averages
avg_temp = calculate_regional_average("data/temperature_file.tif", NILE_DELTA_BBOX)
print(f"Average temperature: {avg_temp:.2f}°C")
```

### Advanced Analysis

```python
# Vapor pressure analysis
analyze_vapor_pressure("data/temperature_file.tif", "results/vapor_pressure.png")

# Regional temperature analysis with visualization
plot_regional_temperature_analysis("data/temperature_file.tif", 
                                  "Nile Delta Temperature Analysis", 
                                  "results/regional_analysis.png")
```

## Data Requirements

The code expects raster data files in GeoTIFF format (.tif) with:
- Temperature data (°C)
- Precipitation data (mm)
- Geographic coordinate system (WGS84 recommended)

### Supported Data Sources

- WorldClim climate data
- Copernicus Climate Change Service (C3S) data
- Other gridded climate datasets in GeoTIFF format

## Configuration

### Geographic Settings

The analysis is configured for the Nile Delta region with these default coordinates:

```python
NILE_DELTA_LON = 31.2  # Longitude
NILE_DELTA_LAT = 30.8  # Latitude

NILE_DELTA_BBOX = {
    "lon_min": 29.5,
    "lon_max": 32.5,
    "lat_min": 30.5,
    "lat_max": 31.5
}
```

You can modify these coordinates for different study areas.

## Dependencies

- `numpy` - Numerical computations
- `matplotlib` - Plotting and visualization
- `rasterio` - Geospatial raster data I/O
- `geopandas` - Geospatial data analysis (optional)

## Output

The code generates:
- High-resolution plots (300 DPI) saved as PNG files
- Statistical summaries printed to console
- Regional average calculations
- Temperature and precipitation maps

## Examples

### Temperature Analysis Example
```python
# Plot temperature map for January
plot_temperature_map("data/temp_jan.tif", 
                    "January Temperature", 
                    "results/temp_jan.png", 
                    "January")
```

### Precipitation Analysis Example
```python
# Analyze precipitation patterns
plot_precipitation_map("data/precip_dec.tif", 
                      "December Precipitation", 
                      "results/precip_dec.png")
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Commit your changes (`git commit -am 'Add new analysis function'`)
4. Push to the branch (`git push origin feature/new-analysis`)
5. Create a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this code in your research, please cite:

```
Sobh, H. (2025). Impact of Climate Change on Rice crop in Kafr Elsheikh region and adaptation measures. 
Master's Thesis, Vrije Universiteit Brussel & Universiteit Gent.
```

## Contact

- **Author**: Hamza Sobh
- **Institution**: Vrije Universiteit Brussel & Universiteit Gent
- **Project**: Master of Science in Sustainable Land Management Engineering

## Acknowledgments

- Prof. Elga Salvadore (Promoter)
- Dr. Ahmed Elnaggar (Co-promotor)
- WorldClim for climate data
- Copernicus Climate Change Service (C3S)

---

*This repository is part of a Master's thesis research on climate change impacts on rice production in Egypt's Nile Delta region.*

