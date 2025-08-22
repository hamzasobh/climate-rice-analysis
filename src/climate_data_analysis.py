#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Climate Data Analysis for Rice Production in Kafr Elsheikh
=========================================================

This script contains functions for analyzing climate data (temperature, precipitation)
and their impacts on rice production in the Kafr Elsheikh region of Egypt.

Author: Hamza Sobh
Date: 2025
Project: Impact of Climate Change on Rice crop in Kafr Elsheikh region
"""

import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.transform import rowcol
from rasterio.mask import mask
import geopandas as gpd

# Geographic coordinates for the Nile Delta region
NILE_DELTA_LON = 31.2  # Longitude
NILE_DELTA_LAT = 30.8  # Latitude

# Bounding box for Nile Delta region
NILE_DELTA_BBOX = {
    "lon_min": 29.5,
    "lon_max": 32.5,
    "lat_min": 30.5,
    "lat_max": 31.5
}

def saturated_vapor_pressure(T):
    """
    Calculate the saturated vapor pressure (e_s) using the Magnus-Tetens formula.

    Parameters:
    -----------
    T : float or numpy array
        Temperature in degrees Celsius.

    Returns:
    --------
    float or numpy array
        Saturated vapor pressure in kPa.
    """
    # Constants for the Magnus-Tetens formula
    A = 0.611  # kPa (converted from hPa by dividing by 10)
    B = 17.27
    C = 237.3  # °C

    # Magnus-Tetens equation
    e_s = A * np.exp((B * T) / (C + T))
    return e_s

def extract_point_value(raster_path, lon, lat):
    """
    Extract raster value at specific geographic coordinates.
    
    Parameters:
    -----------
    raster_path : str
        Path to the raster file
    lon : float
        Longitude coordinate
    lat : float
        Latitude coordinate
        
    Returns:
    --------
    float
        Raster value at the specified coordinates
    """
    with rasterio.open(raster_path) as src:
        raster_data = src.read(1)
        transform = src.transform
        
        # Convert geographic coordinates to pixel indices
        row, col = rowcol(transform, lon, lat)
        
        # Extract value at the pixel
        value = raster_data[row, col]
        
    return value

def plot_temperature_map(raster_path, title, output_path=None, month=None):
    """
    Plot temperature raster data for the Nile Delta region.
    
    Parameters:
    -----------
    raster_path : str
        Path to the temperature raster file
    title : str
        Title for the plot
    output_path : str, optional
        Path to save the plot
    month : str, optional
        Month name for the plot
    """
    with rasterio.open(raster_path) as src:
        raster_data = src.read(1)
        transform = src.transform
        
        # Mask invalid values
        raster_data_masked = np.ma.masked_where(raster_data <= -1e10, raster_data)
        
        # Convert geographic coordinates to pixel indices
        nile_row, nile_col = rowcol(transform, NILE_DELTA_LON, NILE_DELTA_LAT)
        
        # Extract temperature value at Nile Delta
        nile_temp = raster_data[nile_row, nile_col]
        
        # Create plot
        plt.figure(figsize=(10, 8))
        plt.imshow(raster_data_masked, cmap="viridis", vmin=0, vmax=30)
        plt.colorbar(label='Temperature (°C)')
        plt.title(f"{title} - {month}" if month else title)
        
        # Add marker for Nile Delta
        plt.scatter(nile_col, nile_row, color='red', label='Nile Delta', s=30)
        plt.text(nile_col, nile_row - 8, f"{nile_temp:.2f}°C", color="black", fontsize=8)
        plt.legend(loc="lower left")
        
        # Set zoom around Egypt
        zoom_size = 300
        plt.xlim(nile_col - zoom_size, nile_col + zoom_size)
        plt.ylim(nile_row + zoom_size, nile_row - zoom_size)
        
        # Save plot if output path provided
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        
        plt.show()

def plot_precipitation_map(raster_path, title, output_path=None):
    """
    Plot precipitation raster data for the Nile Delta region.
    
    Parameters:
    -----------
    raster_path : str
        Path to the precipitation raster file
    title : str
        Title for the plot
    output_path : str, optional
        Path to save the plot
    """
    with rasterio.open(raster_path) as src:
        precip_data = src.read(1)
        transform = src.transform
        
        # Mask invalid values
        precip_data_masked = np.ma.masked_where(precip_data <= -100, precip_data)
        
        # Convert geographic coordinates to pixel indices
        nile_row, nile_col = rowcol(transform, NILE_DELTA_LON, NILE_DELTA_LAT)
        
        # Extract precipitation value
        nile_precip = precip_data[nile_row, nile_col]
        
        # Create plot
        plt.figure(figsize=(10, 8))
        
        # Set zoom around Egypt
        zoom_size = 80
        plt.xlim(nile_col - zoom_size, nile_col + zoom_size)
        plt.ylim(nile_row + zoom_size, nile_row - zoom_size / 2)
        
        plt.imshow(precip_data_masked, cmap="viridis_r", vmin=0, vmax=15)
        plt.colorbar(label="Precipitation (mm)")
        plt.scatter(nile_col, nile_row, color='black', label='Nile Delta', s=30)
        plt.legend(loc="lower left")
        plt.text(nile_col, nile_row - 8, f"{nile_precip:.2f} mm", color="orange", fontsize=6)
        plt.title(title)
        
        # Save plot if output path provided
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        
        plt.show()
        
        print(f"Precipitation at Nile Delta: {nile_precip:.2f} mm")

def calculate_regional_average(raster_path, bbox):
    """
    Calculate average value for a specific region defined by bounding box.
    
    Parameters:
    -----------
    raster_path : str
        Path to the raster file
    bbox : dict
        Bounding box with keys: lon_min, lon_max, lat_min, lat_max
        
    Returns:
    --------
    float
        Average value for the region
    """
    with rasterio.open(raster_path) as src:
        raster_data = src.read(1)
        transform = src.transform
        
        # Replace nodata with NaN
        raster_data = np.where(raster_data == src.nodata, np.nan, raster_data)
        
        # Convert bounding box to pixel coordinates
        col_min, row_min = ~transform * (bbox["lon_min"], bbox["lat_max"])
        col_max, row_max = ~transform * (bbox["lon_max"], bbox["lat_min"])
        
        # Ensure indices are integers
        row_min, row_max = int(row_min), int(row_max)
        col_min, col_max = int(col_min), int(col_max)
        
        # Crop raster for the region
        region_data = raster_data[row_min:row_max, col_min:col_max]
        
        # Calculate average
        avg_value = np.nanmean(region_data)
        
    return avg_value

def plot_regional_temperature_analysis(raster_path, title, output_path=None):
    """
    Plot temperature analysis for the Nile Delta region with regional average.
    
    Parameters:
    -----------
    raster_path : str
        Path to the temperature raster file
    title : str
        Title for the plot
    output_path : str, optional
        Path to save the plot
    """
    with rasterio.open(raster_path) as src:
        raster_data = src.read(1)
        transform = src.transform
        
        # Replace nodata with NaN
        raster_data = np.where(raster_data == src.nodata, np.nan, raster_data)
        
        # Define display bounding box (larger for visualization)
        display_bbox = {
            "lon_min": 28.0,
            "lon_max": 36.0,
            "lat_min": 26.0,
            "lat_max": 33.0
        }
        
        # Convert bounding boxes to pixel coordinates
        col_min_disp, row_min_disp = ~transform * (display_bbox["lon_min"], display_bbox["lat_max"])
        col_max_disp, row_max_disp = ~transform * (display_bbox["lon_max"], display_bbox["lat_min"])
        
        row_min_disp, row_max_disp = int(row_min_disp), int(row_max_disp)
        col_min_disp, col_max_disp = int(col_min_disp), int(col_max_disp)
        
        # Crop raster for display
        display_data = raster_data[row_min_disp:row_max_disp, col_min_disp:col_max_disp]
        
        # Calculate regional average for Nile Delta
        avg_temp = calculate_regional_average(raster_path, NILE_DELTA_BBOX)
        
        # Create grids for visualization
        cols, rows = np.meshgrid(np.arange(col_min_disp, col_max_disp), 
                                np.arange(row_min_disp, row_max_disp))
        lons, lats = rasterio.transform.xy(transform, rows, cols, offset='center')
        
        # Convert to 2D grids
        lons = np.array(lons).reshape(display_data.shape)
        lats = np.array(lats).reshape(display_data.shape)
        
        # Create plot
        plt.figure(figsize=(10, 8))
        plt.pcolormesh(lons, lats, display_data, cmap="viridis", 
                      shading='auto', vmin=10, vmax=35)
        plt.colorbar(label="Temperature (°C)")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        
        # Set display limits
        plt.xlim(28, 36)
        plt.ylim(26, 33)
        
        # Add average temperature text
        plt.figtext(0.5, -0.05, 
                   f"Average Temperature for the Nile Delta Region: {avg_temp:.2f} °C", 
                   wrap=True, horizontalalignment='center', fontsize=12, style='italic')
        
        plt.title(title)
        plt.tight_layout()
        
        # Save plot if output path provided
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        plt.show()

def analyze_vapor_pressure(temp_raster_path, output_path=None):
    """
    Analyze water vapor pressure at the Nile Delta location.
    
    Parameters:
    -----------
    temp_raster_path : str
        Path to the temperature raster file
    output_path : str, optional
        Path to save the plot
    """
    with rasterio.open(temp_raster_path) as src:
        temp_data = src.read(1)
        transform = src.transform
        
        # Mask invalid values
        temp_data_masked = np.ma.masked_where(temp_data < -100, temp_data)
        
        # Convert coordinates to pixel indices
        nile_row, nile_col = rowcol(transform, NILE_DELTA_LON, NILE_DELTA_LAT)
        
        # Extract temperature at Nile Delta
        nile_temp = temp_data[nile_row, nile_col]
        
        # Calculate vapor pressure
        nile_vapor_pressure = saturated_vapor_pressure(nile_temp)
        
        # Create plot
        plt.figure(figsize=(10, 8))
        
        # Set zoom around region
        zoom_size = 100
        plt.xlim(nile_col - zoom_size, nile_col + zoom_size)
        plt.ylim(nile_row + zoom_size, nile_row - zoom_size)
        
        plt.imshow(temp_data_masked, cmap="plasma_r", vmin=-1, vmax=1)
        plt.colorbar(label="WVP (kPa)")
        plt.scatter(nile_col, nile_row, color='red', label='Nile Delta', s=30)
        plt.legend(loc="lower left")
        plt.text(nile_col, nile_row - 8, f"{nile_vapor_pressure:.3f} kPa", 
                color="white", fontsize=8)
        plt.title("Vapor Pressure at Nile Delta")
        
        # Save plot if output path provided
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        
        plt.show()
        
        print(f"Temperature at Nile Delta: {nile_temp:.2f} °C")
        print(f"Vapor Pressure at Nile Delta: {nile_vapor_pressure:.3f} kPa")

if __name__ == "__main__":
    # Example usage
    print("Climate Data Analysis for Kafr Elsheikh Rice Production")
    print("=" * 55)
    
    # Example file paths (update these with your actual file paths)
    # temp_file = "path/to/temperature/file.tif"
    # precip_file = "path/to/precipitation/file.tif"
    
    # Uncomment and modify these lines to run analysis:
    # plot_temperature_map(temp_file, "Nile Delta Temperature", "output/temp_map.png", "January")
    # plot_precipitation_map(precip_file, "Nile Delta Precipitation", "output/precip_map.png")
    # analyze_vapor_pressure(temp_file, "output/vapor_pressure.png")
    
    print("Analysis functions ready. Update file paths and uncomment function calls to run.")

