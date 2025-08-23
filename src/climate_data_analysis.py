# -*- coding: utf-8 -*-

# import rasterio

# # Provide the correct file path to the .tif file
# file_path = file_path = "/Users/hamza/Documents/Python_Data_2.0/wc2.1_2.5m_tmin-1970-1979/wc2.1_2.5m_tmin_1970_01.tif"
# # Open the .tif file
# with rasterio.open("/Users/hamza/Documents/Python_Data_2.0/wc2.1_2.5m_tmin-1970-1979/wc2.1_2.5m_tmin_1970_01.tif") as src:
#     print("File Metadata:")
#     print(src.meta)  # Print metadata of the .tif file
#     raster_data = src.read(1)  # Read the first band of the raster

   
# # Precipitation data 
# import rasterio
# from rasterio.transform import rowcol
# import numpy as np
# import matplotlib.pyplot as plt

# # File path to the precipitation raster file
# precip_file_path = file_path = "/Users/hamza/Documents/Python_Data_2.0/wc2.0_30s_tavg/wc2.1_30s_tavg_01.tif"

# # Geographic coordinates for the Nile Delta (approximate)
# nile_delta_lon = 31.2  # Longitude
# nile_delta_lat = 30.8  # Latitude

# # Open the precipitation raster file
# with rasterio.open(precip_file_path) as src:
#     precip_data = src.read(1)  # Read the first band of the raster
#     transform = src.transform   # Get the transform (affine)

#     # Mask invalid or extreme values
#     # Mask invalid or non-positive values
#     precip_data_masked = np.ma.masked_where(precip_data <= -100, precip_data)

    
#     # Convert geographic coordinates to pixel indices
#     nile_row, nile_col = rowcol(transform, nile_delta_lon, nile_delta_lat)
    
#     # Extract the precipitation value at the pixel
#     nile_precip = precip_data[nile_row, nile_col]
    
#     # Set plot limits to zoom in around Egypt
#     zoom_size = 80 # Number of pixels to include around Egypt
#     plt.xlim(nile_col - zoom_size, nile_col + zoom_size)
#     plt.ylim(nile_row + zoom_size, nile_row - zoom_size / 2) # Reverse y-axis for raster orientation

#     # Print the precipitation value
#     print(f"Precipitation at Nile Delta (Lon: {nile_delta_lon}, Lat: {nile_delta_lat}): {nile_precip:.2f} mm")

#     # Plot the precipitation raster
#     plt.imshow(precip_data_masked, cmap="viridis_r", vmin=0, vmax=15)  # Use a suitable color map for precipitation
#     plt.colorbar(label="Precipitation (mm)")
#     plt.scatter(nile_col, nile_row, color='black', label='Nile Delta', s=30)
#     plt.legend(loc="lower left")
#     plt.text(nile_col + 0, nile_row - 8, f"{nile_precip:.2f} mm", color="orange", fontsize=6)
#     plt.title("Precip at Nile Delta 2050 - Dec- RCP4.5")
#     # Save the zoomed plot as an image
#     output_plot_path = r"C:\Users\VUB IRMO\Desktop\CC-results\RCP45-Prec-2050\Delta_Prec-Dec2050-RCP4.5.png"
#     plt.savefig(output_plot_path, dpi=300, bbox_inches="tight")  # Save with high resolution
#     plt.show()


# Egypt Map
# import numpy as np
# import matplotlib.pyplot as plt
# import rasterio
# from rasterio.transform import rowcol


# # File path to the raster
# file_path ="/Users/hamza/Documents/Python_Data_2.0/wc2.0_2.5m_tavg/wc2.1_2.5m_tavg_02.tif"

# # Open the raster file and read metadata
# with rasterio.open(file_path) as src:
#     raster_data = src.read(1)  # Read the first band of the raster
#     transform = src.transform  # Get the transform (affine)
#     print("Temp Graph:", raster_data.shape)
#     print("Max Temp:", raster_data.max())
#     print("Min Temp:", raster_data.min())
# # # Replace nodata values with NaN
# # nodata_value = -3.4e+38  # Nodata value from metadata
# # raster_data = raster_data.astype(np.float32)  # Ensure data is in float32 format
# # raster_data[raster_data <= nodata_value] = np.nan  # Replace all nodata values with NaN

# # # Mask invalid values for plotting
# # raster_data_masked = np.ma.masked_invalid(raster_data)

# # Mask invalid values (e.g., very large or very small values)
# raster_data_masked = np.ma.masked_where(raster_data <= -1e10, raster_data)

# # Geographic coordinates of Egypt (approximately)
# egypt_lon = 31.2  # Longitude
# egypt_lat = 30.8  # Latitude



# # Calculate pixel indices for Egypt's location
# egypt_row, egypt_col = rowcol(transform, egypt_lon, egypt_lat)
# # Extract the temperature value at the pixel
# nile_temp = raster_data[egypt_row, egypt_col]
# # Display the masked raster data
# plt.imshow(raster_data_masked, cmap="viridis", vmin=0, vmax=30)
# plt.colorbar(label='Temperature (°C)')
# plt.title("Nile DeltaTemp Feb/1970 - 2000")

# # Add a marker for Egypt using calculated pixel indices
# plt.scatter(egypt_col, egypt_row, color='red', label='Nile Delta', s=30)  # Red marker
# plt.text(egypt_col + 0, egypt_row - 8, f"{nile_temp:.2f}°C", color="black", fontsize=8)
# plt.legend(loc="lower left")
# # plt.pcolormesh(egypt_lon, egypt_lat, raster_data, cmap="viridis", shading='auto', vmin=10, vmax=30)


# # Set plot limits to zoom in around Egypt
# zoom_size = 300  # Number of pixels to include around Egypt
# plt.xlim(egypt_col - zoom_size, egypt_col + zoom_size)
# plt.ylim(egypt_row + zoom_size, egypt_row - zoom_size)  # Reverse y-axis for raster orientation

# # Save the zoomed plot as an image
# output_plot_path = file_path = "/Users/hamza/Documents/cc-results/NileDeltaTempFeb.tif"
# plt.savefig(output_plot_path, dpi=300, bbox_inches="tight")  # Save with high resolution

# # Show the zoomed plot
# plt.show()

# # The Nile Delta region 
# import numpy as np
# import matplotlib.pyplot as plt
# from rasterio.transform import rowcol
# import rasterio

# # File path to the raster file
# file_path = r"C:\Users\VUB IRMO\Desktop\Python\wc2.1_10m_tmin\wc2.1_10m_tmin_01.tif"

# # Open the raster file
# with rasterio.open(file_path) as src:
#     raster_data = src.read(1)  # Read the first band
#     transform = src.transform  # Get the transform for row/col conversion

# # Mask invalid values
# raster_data_masked = np.ma.masked_where(raster_data <= -1e10, raster_data)

# # Geographic coordinates for the Nile Delta region
# lon_min, lon_max = 29.0, 33.0  # Longitudes
# lat_min, lat_max = 30.0, 32.0  # Latitudes

# # Calculate pixel coordinates from geographic coordinates
# row_min, col_min = rowcol(transform, lon_min, lat_max)  # Top-left
# row_max, col_max = rowcol(transform, lon_max, lat_min)  # Bottom-right

# # Display the zoomed-in raster data
# plt.imshow(raster_data_masked[row_min:row_max, col_min:col_max], 
#             cmap="viridis", extent=[lon_min, lon_max, lat_min, lat_max])
# plt.colorbar(label="Temperature (°C)")
# plt.title(" Nile Delta Region (Jan)")

# # Save the zoomed plot as an image
# # output_plot_path = r"C:\Users\VUB IRMO\Desktop\Python\Temp_Hist\Delta_Temp-Jan_T= 9.315999984741211°C.png"
# # plt.savefig(output_plot_path, dpi=300, bbox_inches="tight")  # Save with high resolution

# # Add grid lines for clarity
# plt.grid(True)

# # Show the plot
# plt.show()


# # Egypt Temperature
# import rasterio
# from rasterio.transform import rowcol

# # File path to the raster
# file_path = r"C:\Users\VUB IRMO\Desktop\Python\wc2.1_10m_tmin\wc2.1_10m_tmin_12.tif"

# # Geographic coordinates of Egypt (approximate)
# egypt_lon = 30.0  # Longitude
# egypt_lat = 26.0  # Latitude

# # Open the raster file
# with rasterio.open(file_path) as src:
#     raster_data = src.read(1)  # Read the first band of the raster
#     transform = src.transform  # Get the transform (affine)
    
#     # Convert geographic coordinates to pixel indices
#     egypt_row, egypt_col = rowcol(transform, egypt_lon, egypt_lat)
    
#     # Extract the temperature value at the pixel
#     egypt_temp = raster_data[egypt_row, egypt_col]

#     # Print the temperature value
#     print(f"Temperature at Egypt (Lon: {egypt_lon}, Lat: {egypt_lat}): {egypt_temp}°C")


# The Nile Delta Temperature 

# import rasterio
# from rasterio.transform import rowcol

# # File path to the raster
# file_path = r"C:\Users\VUB IRMO\Desktop\Python\wc2.1_10m_tmin\wc2.1_10m_tmin_12.tif"

# # Geographic coordinates for a specific location in the Nile Delta (example coordinates)
# nile_delta_lon = 31.2  # Longitude (central Nile Delta)
# nile_delta_lat = 30.8  # Latitude (central Nile Delta)

# # Open the raster file
# with rasterio.open(file_path) as src:
#     raster_data = src.read(1)  # Read the first band of the raster
#     transform = src.transform  # Get the transform (affine)

#     # Convert geographic coordinates to pixel indices
#     nile_row, nile_col = rowcol(transform, nile_delta_lon, nile_delta_lat)

#     # Extract the temperature value at the pixel
#     nile_temp = raster_data[nile_row, nile_col]

#     # Print the temperature value
#     print(f"Temperature at Nile Delta (Lon: {nile_delta_lon}, Lat: {nile_delta_lat}): {nile_temp:.2f}°C")







"""
from rasterio import Affine
from rasterio.enums import Resampling

# Define output file path
output_path = r"C:/Users/VUB IRMO/Desktop/Python/temp.tif"

# Save modified raster
with rasterio.open(
    output_path,
    'w',
    driver='GTiff',
    height=raster_data.shape[0],
    width=raster_data.shape[1],
    count=1,  # Number of bands
    dtype=raster_data.dtype,
    crs=src.crs,  # Coordinate Reference System
    transform=src.transform,  # Geo-transform
) as dst:
    dst.write(raster_data, 1)
    """
    # /Users/hamza/Documents/Python_Data_2.0/wc2.0_2.5m_tavg/wc2.1_2.5m_tavg_02.tif
# import numpy as np
# import rasterio
# import matplotlib.pyplot as plt
# import geopandas as gpd
# from rasterio.plot import show

# # File path to the raster
# file_path = "/Users/hamza/Documents/Python_Data_2.0/wc2.1_2.5m_tmax-1970-1979/wc2.1_2.5m_tmax_1970_12.tif"

# # Open the raster file
# with rasterio.open(file_path) as src:
#     raster_data = src.read(1)
#     transform = src.transform
#     raster_data = np.where(raster_data == src.nodata, np.nan, raster_data)
    
    

#     # Define bounding box for Nile Delta
    
#     bounding_box_temp = {  
#     "lon_min": 29.5,  # Keep original values
#     "lon_max": 32.5,  
#     "lat_min": 30.5,  
#     "lat_max": 31.5  
# }
    
#     bounding_box_display = {  
#     "lon_min": 28.0,  # Expanded for visualization
#     "lon_max": 36.0,  
#     "lat_min": 26.0,  
#     "lat_max": 33.0  
# }

#     # Convert bounding box to pixel coordinates
#     col_min_temp, row_min_temp = ~transform * (bounding_box_temp["lon_min"], bounding_box_temp["lat_max"])
#     col_max_temp, row_max_temp = ~transform * (bounding_box_temp["lon_max"], bounding_box_temp["lat_min"])
    
#     col_min_disp, row_min_disp = ~transform * (bounding_box_display["lon_min"], bounding_box_display["lat_max"])
#     col_max_disp, row_max_disp = ~transform * (bounding_box_display["lon_max"], bounding_box_display["lat_min"])

#       # Ensure indices are integers
#     row_min_temp, row_max_temp = int(row_min_temp), int(row_max_temp)
#     col_min_temp, col_max_temp = int(col_min_temp), int(col_max_temp)

#     row_min_disp, row_max_disp = int(row_min_disp), int(row_max_disp)
#     col_min_disp, col_max_disp = int(col_min_disp), int(col_max_disp)

#     # Crop raster for temperature calculation
#     temp_data = raster_data[row_min_temp:row_max_temp, col_min_temp:col_max_temp]

#     # Crop raster for display (larger region)
#     display_data = raster_data[row_min_disp:row_max_disp, col_min_disp:col_max_disp]

#     # Calculate average temperature for original box
#     avg_temp = np.nanmean(temp_data)

#     # Create grids for visualization
#     cols, rows = np.meshgrid(np.arange(col_min_disp, col_max_disp), np.arange(row_min_disp, row_max_disp))
#     lons, lats = rasterio.transform.xy(transform, rows, cols, offset='center')

# # Convert longitude and latitude arrays into 2D grids
# lons = np.array(lons).reshape(display_data.shape)
# lats = np.array(lats).reshape(display_data.shape)

# # Plot the expanded region but keep original temp calculations
# plt.figure(figsize=(10, 8))
# plt.pcolormesh(lons, lats, display_data, cmap="viridis", shading='auto', vmin=10, vmax=35)
# plt.colorbar(label="Temperature (°C)")
# plt.xlabel("Longitude")
# plt.ylabel("Latitude")

# # Set display limits to show the expanded area
# plt.xlim(28, 36)  # Show the larger longitude range
# plt.ylim(26, 33)  # Show the larger latitude range



# # Add the average temperature text below the plot
# plt.figtext(0.5, -0.05,  # Adjusted position
#             f"Maximum Temperature for the Nile Delta Region: {avg_temp:.2f} °C", 
#             wrap=True, horizontalalignment='center', fontsize=12, style='italic')

# # Add a title
# plt.title("Nile Delta Temperature - Dec 1970 ")

# # Save the plot to a file
# output_file = "/Users/hamza/Documents/CC-results/HistoricalTemp/maxTemp-1970/NileDeltaTemperature-Dec.png"
# plt.tight_layout()  # Adjust layout to fit everything
# plt.savefig(output_file, dpi=300, bbox_inches='tight')  # Save with high resolution


# plt.show()








# # Final code comparison chart
# import numpy as np
# import rasterio
# import geopandas as gpd
# import matplotlib.pyplot as plt
# from rasterio.mask import mask

# # Define file paths for min, max, and avg temperature rasters
# file_path_min_template = "/Users/hamza/Documents/Python_Data_2.1/wc2.1_2.5m_tmin/wc2.1_2.5m_tmin_{month}.tif"
# file_path_max_template = "/Users/hamza/Documents/Python_Data_2.1/wc2.1_2.5m_tmax/wc2.1_2.5m_tmax_{month}.tif"
# file_path_avg_template = "/Users/hamza/Documents/Python_Data_2.1/wc2.1_2.5m_tavg/wc2.1_2.5m_tavg_{month}.tif"

# # Shapefile of the Nile Delta region
# shapefile_path = "/Users/hamza/Documents/Egypt Data/EGY_adm/EGY_adm2.shp"



# # Load the shapefile
# gdf = gpd.read_file(shapefile_path)

# # Define the list of Nile Delta governorates
# nile_delta_governorates = [
#     "Al Buhayrah", "Kafr ash Shaykh", "Ad Daqahliyah",
#     "Al Gharbiyah", "Al Minufiyah", "Ash Sharqiyah", "Dumyat"
# ]

# # Filter the shapefile to keep only the Nile Delta regions
# gdf_nile_delta = gdf[gdf["NAME_1"].isin(nile_delta_governorates)]



# # Filter for Nile Delta (Ensure your shapefile has a column that identifies the region)
# gdf_nile_delta = gdf[gdf["NAME_1"].str.contains("Nile Delta", case=False, na=False)]  # Adjust column name


# # Function to extract temperature data using shapefile mask
# def extract_temperature(file_path, gdf_nile_delta):
#     try:
#         with rasterio.open(file_path) as src:
#             out_image, out_transform = mask(src, gdf_nile_delta.geometry, crop=True)
#             raster_data = np.where(out_image == src.nodata, np.nan, out_image)

#             # Compute mean temperature for the region (ignoring NaN values)
#             return np.nanmean(raster_data)

#     except Exception as e:
#         print(f"Error processing file {file_path}: {e}")
#         return np.nan  # Return NaN if an error occurs

# # Initialize lists to store temperature data
# min_temps = []
# max_temps = []
# avg_temps = []

# months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# # Loop through each month, extract temperature values
# for month in range(1, 13):
#     file_path_min = file_path_min_template.format(month=f"{month:02d}")
#     file_path_max = file_path_max_template.format(month=f"{month:02d}")
#     file_path_avg = file_path_avg_template.format(month=f"{month:02d}")

#     min_temp = extract_temperature(file_path_min, gdf)
#     max_temp = extract_temperature(file_path_max, gdf)
#     avg_temp = extract_temperature(file_path_avg, gdf)

#     min_temps.append(min_temp)
#     max_temps.append(max_temp)
#     avg_temps.append(avg_temp)

# # Plot the temperature comparison chart
# plt.figure(figsize=(10, 6))
# plt.plot(months, min_temps, marker='o', label="Min Temperature", color="blue")
# plt.plot(months, max_temps, marker='o', label="Max Temperature", color="red")
# plt.plot(months, avg_temps, marker='o', label="Avg Temperature", color="green")

# # Add temperature values on the chart
# for i, month in enumerate(months):
#     plt.text(month, min_temps[i] + 0.7, f"{min_temps[i]:.1f}°C", ha='center', fontsize=6, color="blue", weight="bold", bbox=dict(facecolor='white', alpha=0.6, edgecolor='blue', boxstyle='round,pad=0.2'))
#     plt.text(month, max_temps[i] + 0.8, f"{max_temps[i]:.1f}°C", ha='center', fontsize=6, color="red", weight="bold", bbox=dict(facecolor='white', alpha=0.6, edgecolor='blue', boxstyle='round,pad=0.2'))
#     plt.text(month, avg_temps[i] + 0.7, f"{avg_temps[i]:.1f}°C", ha='center', fontsize=6, color="green", weight="bold", bbox=dict(facecolor='white', alpha=0.6, edgecolor='blue', boxstyle='round,pad=0.2'))


# plt.xlabel("Months")
# plt.ylabel("Temperature (°C)")
# plt.title("Comparison of Min, Max, and Avg Temperature for Nile Delta (1970-2000)")
# plt.legend(loc="upper right", fontsize=8, markerscale=0.8, frameon=True, edgecolor="black")
# plt.grid(False)
# # Define output file path
# output_chart = "/Users/hamza/Documents/CC-results/HistoricalTemp/Comparison_Charts/Comparison_Temperature_Nile_Delta.png"

# # Save the figure
# plt.savefig(output_chart, dpi=300, bbox_inches='tight')

# # Print confirmation message
# print(f"Temperature comparison chart saved to: {output_chart}")

# plt.show()




##Average_Temp
# import numpy as np
# import rasterio
# import geopandas as gpd
# import matplotlib.pyplot as plt
# from rasterio.mask import mask
# from rasterio.plot import show

# # File paths (Update these paths)
# raster_path = "/Users/hamza/Documents/Python_Data_2.0/wc2.0_2.5m_tavg/wc2.1_2.5m_tavg_02.tif"
# shapefile_path = "/Users/hamza/Documents/Egypt Data/EGY_adm/EGY_adm2.shp"

# # Load the shapefile (Nile Delta region)
# gdf = gpd.read_file(shapefile_path)

# # Open the raster file
# with rasterio.open(raster_path) as src:
#     raster_data = src.read(1)  # Read the first band
#     raster_crs = src.crs  # Get CRS
#     transform = src.transform  # Get transformation



#     # Ensure the shapefile is in the same CRS as the raster
#     if gdf.crs != raster_crs:
#         gdf = gdf.to_crs(raster_crs)
        

#     # Mask raster using the shapefile (extract only within boundaries)
#     out_image, out_transform = mask(src, gdf.geometry, crop=True)
    
#     # Replace no-data values with NaN
#     out_image = np.where(out_image == src.nodata, np.nan, out_image)

# # Compute statistics for the Nile Delta region
# mean_temp = np.nanmean(out_image)
# min_temp = np.nanmin(out_image)
# max_temp = np.nanmax(out_image)

# # Print temperature statistics
# print(f"Nile Delta Temperature Statistics:")
# print(f"  Mean Temperature: {mean_temp:.2f}°C")
# print(f"  Min Temperature: {min_temp:.2f}°C")
# print(f"  Max Temperature: {max_temp:.2f}°C")


# # Plot the extracted temperature region
# fig, ax = plt.subplots(figsize=(10, 6))
# show(out_image, transform=out_transform, ax=ax, cmap="coolwarm", alpha=0.7)
# gdf.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1.5)
# ax.set_title("Extracted Temperature for Nile Delta Region", fontsize=14)
# plt.show()




# #TemperatureNileDeltaMap
# import numpy as np
# import geopandas as gpd
# import matplotlib.pyplot as plt
# import rasterio
# from rasterio.mask import mask
# from rasterio.plot import plotting_extent
# import os

# # === FILE PATHS ===
# raster_path = "/Users/hamza/Documents/CC-results/Future-Data/wc2.1_2.5m_tmin_ACCESS-CM2_ssp585_2021-2040.tif"  # 🌡️ Replace with temperature raster
# shapefile_path = "/Users/hamza/Documents/Thesis/Egypt Data/EGY_adm/EGY_adm1.shp"
# output_file = "/Users/hamza/Documents/CC-results/HistoricalTemp/Futuremax.png"

# # === Ensure output folder exists ===
# os.makedirs(os.path.dirname(output_file), exist_ok=True)

# # === Load shapefile ===
# gdf = gpd.read_file(shapefile_path)

# # === Filter governorates ===
# nile_delta_governorates = [
#     "Al Buhayrah", "Kafr ash Shaykh", "Ad Daqahliyah",
#     "Al Gharbiyah", "Al Minufiyah", "Ash Sharqiyah", "Dumyat"
# ]
# gdf_display = gdf[gdf["NAME_1"].isin(nile_delta_governorates)]
# gdf_kafr = gdf[gdf["NAME_1"] == "Kafr ash Shaykh"]

# # === Open raster and process ===
# with rasterio.open(raster_path) as src:
#     if gdf_display.crs != src.crs:
#         gdf_display = gdf_display.to_crs(src.crs)
#         gdf_kafr = gdf_kafr.to_crs(src.crs)

#     # === Clip raster to Kafr ash Shaykh (for value extraction only) ===
#     clipped_image, _ = mask(src, gdf_kafr.geometry, crop=True)
#     nodata_val = src.nodata
#     temp_data = np.where(clipped_image == nodata_val, np.nan, clipped_image)
#     valid_vals = temp_data[np.isfinite(temp_data)]

#     if valid_vals.size > 0:
#         mean_temp = np.nanmean(valid_vals)
#         print(f"📊 Mean Temperature in Kafr ash Shaykh: {mean_temp:.2f} °C")

#     else:
#         print("⚠️ No valid temperature values found in Kafr ash Shaykh.")

#     # === Clip raster to Nile Delta for display ===
#     display_image, display_transform = mask(src, gdf_display.geometry, crop=True)
#     display_array = np.where(display_image == nodata_val, np.nan, display_image)
#     extent = plotting_extent(display_array[0], display_transform)

# # === Plot ===
# fig, ax = plt.subplots(figsize=(12, 8))
# im = ax.imshow(display_array[0], extent=extent, cmap='OrRd', vmin=0, vmax=35)

# # === Boundaries & Labels ===
# gdf_display.boundary.plot(ax=ax, edgecolor='black', linewidth=1.5)
# for _, row in gdf_display.iterrows():
#     centroid = row.geometry.centroid
#     ax.text(centroid.x, centroid.y, row["NAME_1"], fontsize=7, ha='center', color='black')

# # === Colorbar ===
# cbar = plt.colorbar(im, ax=ax, orientation="vertical", shrink=0.8, pad=0.03)
# cbar.set_label("Avg Temperature (°C)", fontsize=12)

# # === Title ===
# ax.set_title("Kafr Elsheikh – Avg Temperature (Dec 1970–2000)", fontsize=13)

# ax.text(extent[0] + 0.05, extent[3] - 0.08,
#         f"Avg Temp – Kafr El-Sheikh: {mean_temp:.2f} °C",
#         fontsize=8, color='black',
#         bbox=dict(facecolor='white', edgecolor='black'))

# # === Axes & Grid ===
# ax.set_xlabel("Longitude")
# ax.set_ylabel("Latitude")

# # === Save and show ===
# plt.tight_layout()
# plt.savefig(output_file, dpi=300)
# plt.show()



# #FinalPrecepPlot
# # === Nile Delta Precipitation Map – Jan 1970-2000 ===

# import numpy as np
# import geopandas as gpd
# import matplotlib.pyplot as plt
# import rasterio
# from rasterio.mask import mask
# from rasterio.plot import show, plotting_extent
# import os

# # === FILE PATHS ===
# raster_path = "/Users/hamza/Documents/Python_Data_2.1/wc2.1_2.5m_prec/wc2.1_2.5m_prec_01.tif"
# shapefile_path = "/Users/hamza/Documents/Thesis/Egypt Data/EGY_adm/EGY_adm1.shp"
# output_file = "/Users/hamza/Documents/CC-results/HistoricalPrecip/NileDelta_Precip_Jan_1970_2000_Map.png"

# # === Ensure output folder exists ===
# os.makedirs(os.path.dirname(output_file), exist_ok=True)

# # === Load shapefile ===
# gdf = gpd.read_file(shapefile_path)

# # === Filter only Nile Delta governorates ===
# nile_delta_governorates = [
#     "Al Buhayrah", "Kafr ash Shaykh", "Ad Daqahliyah",
#     "Al Gharbiyah", "Al Minufiyah", "Ash Sharqiyah", "Dumyat"
# ]
# gdf_filtered = gdf[gdf["NAME_1"].isin(nile_delta_governorates)]

# # === Open raster and clip ===
# with rasterio.open(raster_path) as src:
#     if gdf_filtered.crs != src.crs:
#         gdf_filtered = gdf_filtered.to_crs(src.crs)

#     out_image, out_transform = mask(src, gdf_filtered.geometry, crop=True)
#     nodata_val = src.nodata
#     out_image = np.where(out_image == nodata_val, np.nan, out_image)

# # === Extract precipitation values inside polygon ===
# precip_values = out_image[0]
# valid_vals = precip_values[np.isfinite(precip_values)]

# if valid_vals.size > 0:
#     avg_precip = np.mean(valid_vals)
#     print(f"🌧️ Average Jan Precipitation in Nile Delta: {avg_precip:.2f} mm")
# else:
#     print("⚠️ No valid precipitation values found.")
#     avg_precip = 0

# # === Plot ===
# fig, ax = plt.subplots(figsize=(12, 8))
# extent = plotting_extent(precip_values, out_transform)

# im = ax.imshow(precip_values, extent=extent, cmap='YlGnBu', vmin=0, vmax=50)
# gdf_filtered.boundary.plot(ax=ax, edgecolor='black', linewidth=1.5)

# # === Annotate governorates ===
# for idx, row in gdf_filtered.iterrows():
#     centroid = row.geometry.centroid
#     ax.text(centroid.x, centroid.y, row["NAME_1"], fontsize=7, ha='center', color='black')

# # === Colorbar ===
# cbar = plt.colorbar(im, ax=ax, orientation="vertical", shrink=0.8, pad=0.03)
# cbar.set_label("Precipitation (mm)", fontsize=12)

# # === Title + Annotation ===
# ax.set_title("Nile Delta – Avg Precipitation (Jan 1970–2000)", fontsize=14)
# ax.text(extent[0]+0.1, extent[3]-0.3, f"Avg Jan Precip: {avg_precip:.2f} mm",
#         fontsize=10, color='black', bbox=dict(facecolor='white', edgecolor='black'))

# # === Display coordinate grid ===
# ax.set_xlabel("Longitude")
# ax.set_ylabel("Latitude")
# # ax.grid(True, linestyle="--", alpha=0.4)

# # === Save and show ===
# plt.tight_layout()
# plt.savefig(output_file, dpi=300)
# plt.show()



# #Kafr_elsheikh_Prec
# import numpy as np
# import geopandas as gpd
# import matplotlib.pyplot as plt
# import rasterio
# from rasterio.mask import mask
# from rasterio.plot import plotting_extent
# import os

# # === FILE PATHS ===
# raster_path = "/Users/hamza/Documents/Python_Data_2.1/wc2.1_2.5m_prec/wc2.1_2.5m_prec_12.tif"
# shapefile_path = "/Users/hamza/Documents/Thesis/Egypt Data/EGY_adm/EGY_adm1.shp"
# output_file = "/Users/hamza/Documents/CC-results/HistoricalPrecip/Kafr_Elsheikh_Prec_Dec.png"

# # === Ensure output folder exists ===
# os.makedirs(os.path.dirname(output_file), exist_ok=True)

# # === Load shapefile ===
# gdf = gpd.read_file(shapefile_path)

# # === Filter governorates ===
# nile_delta_governorates = [
#     "Al Buhayrah", "Kafr ash Shaykh", "Ad Daqahliyah",
#     "Al Gharbiyah", "Al Minufiyah", "Ash Sharqiyah", "Dumyat"
# ]
# gdf_display = gdf[gdf["NAME_1"].isin(nile_delta_governorates)]
# gdf_kafr = gdf[gdf["NAME_1"] == "Kafr ash Shaykh"]

# # === Open raster and process ===
# with rasterio.open(raster_path) as src:
#     if gdf_display.crs != src.crs:
#         gdf_display = gdf_display.to_crs(src.crs)
#         gdf_kafr = gdf_kafr.to_crs(src.crs)

#     # === Clip raster to Kafr ash Shaykh for stats ===
#     clipped_image, _ = mask(src, gdf_kafr.geometry, crop=True)
#     nodata_val = src.nodata
#     precip_data = np.where(clipped_image == nodata_val, np.nan, clipped_image)
#     valid_vals = precip_data[np.isfinite(precip_data)]

#     if valid_vals.size > 0:
#         min_precip = np.nanmin(valid_vals)
#         max_precip = np.nanmax(valid_vals)
#         mean_precip = np.nanmean(valid_vals)
#         print("📊 Precipitation stats for Kafr ash Shaykh:")
#         print(f"▪️ Minimum: {min_precip:.2f} mm")
#         print(f"▪️ Maximum: {max_precip:.2f} mm")
#         print(f"▪️ Mean   : {mean_precip:.2f} mm")
#     else:
#         print("⚠️ No valid precipitation values found in Kafr ash Shaykh.")
#         min_precip = max_precip = mean_precip = 0

#     # === Clip raster to full Nile Delta for display ===
#     display_image, display_transform = mask(src, gdf_display.geometry, crop=True)
#     display_array = np.where(display_image == nodata_val, np.nan, display_image)
#     extent = plotting_extent(display_array[0], display_transform)

# # === Plot ===
# fig, ax = plt.subplots(figsize=(12, 8))
# im = ax.imshow(display_array[0], extent=extent, cmap='YlGnBu', vmin=0, vmax=55)

# # === Boundaries & Labels ===
# gdf_display.boundary.plot(ax=ax, edgecolor='black', linewidth=1.5)
# for _, row in gdf_display.iterrows():
#     centroid = row.geometry.centroid
#     ax.text(centroid.x, centroid.y, row["NAME_1"], fontsize=7, ha='center', color='black')

# # === Colorbar ===
# cbar = plt.colorbar(im, ax=ax, orientation="vertical", shrink=0.8, pad=0.03)
# cbar.set_label("Precipitation (mm)", fontsize=12)

# # === Title & Annotated Stats ===
# ax.set_title("Kafr Elsheikh – Precipitation (Dec 1970–2000)", fontsize=13)
# ax.text(extent[0] + 0.1, extent[3] - 0.2,
#         f"Kafr El-Sheikh:\nMin: {min_precip:.1f} mm\nMax: {max_precip:.1f} mm\nMean: {mean_precip:.1f} mm",
#         fontsize=8, color='black', bbox=dict(facecolor='white', edgecolor='black'))

# # === Axes & Grid ===
# ax.set_xlabel("Longitude")
# ax.set_ylabel("Latitude")

# # === Save and show ===
# plt.tight_layout()
# plt.savefig(output_file, dpi=300)
# plt.show()



# #TemperatureComparisonFig
# import os
# import numpy as np
# import rasterio
# import matplotlib.pyplot as plt
# import geopandas as gpd
# from rasterio.mask import mask

# # Define file paths
# raster_path_template = "/Users/hamza/Documents/Python_Data_2.1/wc2.1_2.5m_tmin/wc2.1_2.5m_tmin_{month:02d}.tif"  # Change to correct raster type
# shapefile_path = "/Users/hamza/Documents/Thesis/Egypt Data/EGY_adm/EGY_adm1.shp"

# # Load the shapefile
# gdf = gpd.read_file(shapefile_path)



# # Function to extract temperature values using the shapefile mask
# def extract_temperature(raster_path, shapefile):
#     try:
#         with rasterio.open(raster_path) as src:
#             out_image, out_transform = mask(src, shapefile.geometry, crop=True)
#             raster_data = np.where(out_image == src.nodata, np.nan, out_image)
#             return raster_data.flatten()  # Return all extracted values
#     except Exception as e:
#         print(f"Error processing file {raster_path}: {e}")
#         return np.nan  # Return NaN if error occurs

# # Initialize lists for monthly storage
# months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
# real_temps = []
# min_temps = []
# max_temps = []

# # Loop through each month and extract temperature values
# for month in range(1, 13):  # Months 1 to 12
#     raster_path = raster_path_template.format(month=month)  # Format month
    
#     if os.path.exists(raster_path):  # Ensure raster file exists
#         month_temps = extract_temperature(raster_path, gdf) # Extract data

#         # Compute real, min, and max temperature for the Nile Delta
#         real_temps.append(np.nanmean(month_temps))  # Mean real temperature from extracted values
#         min_temps.append(np.nanmin(month_temps))    # Minimum extracted temperature
#         max_temps.append(np.nanmax(month_temps))    # Maximum extracted temperature
#     else:
#         print(f"File not found: {raster_path}")
#         real_temps.append(np.nan)
#         min_temps.append(np.nan)
#         max_temps.append(np.nan)

# # Plot Comparison Chart
# plt.figure(figsize=(10, 6))
# plt.plot(months, min_temps, marker="o", label=" Min-Avg Temp", color="blue")
# plt.plot(months, max_temps, marker="o", label=" Max-Avg Temp", color="red")
# plt.plot(months, real_temps, marker="o", label="Avg Temp", color="green")
  
# # Add temperature values on the chart
# for i, month in enumerate(months):
#     plt.text(month, min_temps[i] + 0.7, f"{min_temps[i]:.1f}°C", ha='center', fontsize=6, color="blue", weight="bold", bbox=dict(facecolor='white', alpha=0.6, edgecolor='blue', boxstyle='round,pad=0.2'))
#     plt.text(month, max_temps[i] + 0.8, f"{max_temps[i]:.1f}°C", ha='center', fontsize=6, color="red", weight="bold", bbox=dict(facecolor='white', alpha=0.6, edgecolor='blue', boxstyle='round,pad=0.2'))
#     plt.text(month, real_temps[i] + 0.7, f"{real_temps[i]:.1f}°C", ha='center', fontsize=6, color="green", weight="bold", bbox=dict(facecolor='white', alpha=0.6, edgecolor='blue', boxstyle='round,pad=0.2'))
# plt.xlabel("Months")
# plt.ylabel("Temperature (°C)")
# plt.title("Comparison of Min,Max, and Avg Min Temperature for Kafr Elsheikh (1970-2000)")

# # Add shaded regions with shredded (hatched) effect
# plt.fill_between(months, min_temps, max_temps, color='red', alpha=0.15, label="_nolegend_")  # Remove from legend
# plt.fill_between(months, min_temps, real_temps, color='cyan', alpha=0.15, label="_nolegend_")  # Remove from legend

# plt.legend()
# plt.legend(loc="upper right", fontsize=8, markerscale=0.8, frameon=True, edgecolor="black")
# plt.savefig("/Users/hamza/Documents/CC-results/HistoricalTemp/Comparison_Charts/Comparison_Avg_Temperatures3.png", dpi=300)
# plt.show()



#NileRiver_Shapefile

# import geopandas as gpd
# import matplotlib.pyplot as plt
# import contextily as ctx

# # Load the shapefile
# shapefile_path = "/Users/hamza/Documents/Python_Data_2.1/EGY_wat/EGY_water_areas_dcw.shp"
# water_bodies = gpd.read_file(shapefile_path)

# # Display attribute columns to find the Nile attribute name
# print(water_bodies.columns)
# print(water_bodies.head())

# # Filter for Nile River (assuming 'NAME' contains names of water bodies)
# nile_river = water_bodies[water_bodies['NAME'].str.contains('Nile', case=False, na=False)]

# # Check the filtered Nile data
# print(nile_river)

# import rasterio
# from rasterio.plot import show

# # Load raster
# raster_path = "/Users/hamza/Documents/Python_Data_2.1/wc2.1_2.5m_tmax/wc2.1_2.5m_tmax_12.tif"
# with rasterio.open(raster_path) as src:
#     fig, ax = plt.subplots(figsize=(10, 10))
#     show(src, ax=ax)  # Show raster as background

#     # Plot Nile (make sure CRS matches raster CRS)
#     nile_river = nile_river.to_crs(src.crs)
#     nile_river.plot(ax=ax, linewidth=2, edgecolor='blue')
    
#     # Set axis limits to Nile River bounds with some padding
#     minx, miny, maxx, maxy = nile_river.total_bounds

#     # Add padding (e.g., 1 or 2 degrees for latitude/longitude)
#     padding_x = 4.0  # adjust this value as needed
#     padding_y = 2.0  # adjust this value as needed

#     # Apply padding
#     ax.set_xlim(minx - padding_x, maxx + padding_x)
#     ax.set_ylim(miny - padding_y, maxy + padding_y)

#     plt.show()




# #NDVI_Plot

# import rasterio
# import numpy as np
# import matplotlib.pyplot as plt
# import rasterio.plot

# # ✅ Define the file path (Make sure it's correct)
# file_path = "/Users/hamza/Downloads/NDVI_Nile_Delta.tif"

# # ✅ Load NDVI GeoTIFF file
# with rasterio.open(file_path) as src:
#     ndvi = src.read(1)  # Read the first (only) band
#     crs = src.crs  # Get CRS (Projection)
#     bounds = src.bounds  # Get spatial extent
#     transform = src.transform  # Get transformation for plotting

# # ✅ Print metadata for verification
# print("CRS:", crs)
# print("Bounds:", bounds)
# print("Min NDVI:", np.min(ndvi))
# print("Max NDVI:", np.max(ndvi))

# # ✅ Convert UTM (meters) to Longitude/Latitude
# from rasterio.warp import transform_bounds

# lon_min, lat_min, lon_max, lat_max = transform_bounds(crs, "EPSG:4326", *bounds)
# print("Geographic Bounds:", lon_min, lat_min, lon_max, lat_max)

# # ✅ Create Latitude & Longitude grid
# height, width = ndvi.shape
# x = np.linspace(lon_min, lon_max, width)  # Longitudes
# y = np.linspace(lat_max, lat_min, height)  # Latitudes (inverted)

# # ✅ Plot NDVI Map (Georeferenced)
# plt.figure(figsize=(10, 6))
# plt.imshow(ndvi, extent=[lon_min, lon_max, lat_min, lat_max], cmap="RdYlGn", vmin=-1, vmax=1)
# cbar = plt.colorbar(label="NDVI Value")
# plt.title("NDVI Map - Nile Delta (Georeferenced)")
# plt.xlabel("Longitude")
# plt.ylabel("Latitude")
# plt.show()

# # ✅ Plot NDVI Histogram
# plt.figure(figsize=(8, 5))
# plt.hist(ndvi.flatten(), bins=50, color="green", alpha=0.7)
# plt.xlabel("NDVI Value")
# plt.ylabel("Frequency")
# plt.title("NDVI Value Distribution - Nile Delta")
# plt.grid()
# plt.show()



# #NDVI_Comparison
# import rasterio
# import numpy as np
# import matplotlib.pyplot as plt
# import rasterio.plot
# from rasterio.warp import transform_bounds

# # ✅ Define NDVI file paths
# ndvi_files = {
#     "January 2023": "/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jan2023.tif",
#     "July 2023": "/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jul2023.tif",
#     "January 2024": "/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jan2024.tif",
#     "July 2024": "/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jul2024.tif"
# }

# # ✅ Function to Load & Process NDVI
# def load_ndvi(file_path):
#     with rasterio.open(file_path) as src:
#         ndvi = src.read(1)  # Read NDVI band
#         crs = src.crs  # Coordinate Reference System
#         bounds = src.bounds  # Raster boundaries
#         transform = src.transform  # Transformation matrix

#     # Convert UTM to Latitude/Longitude
#     lon_min, lat_min, lon_max, lat_max = transform_bounds(crs, "EPSG:4326", *bounds)

#     return ndvi, lon_min, lat_min, lon_max, lat_max

# # ✅ Load all NDVI datasets
# ndvi_data = {date: load_ndvi(file) for date, file in ndvi_files.items()}

# # ✅ Step 2: Plot NDVI Over Time with Enhanced Contrast (Side-by-Side Comparison)
# fig, axes = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid for side-by-side comparison
# axes = axes.flatten()  # Convert to 1D array

# # ✅ Set improved contrast range for NDVI
# vmin, vmax = -0.2, 0.6  # Balanced NDVI contrast

# for ax, (date, (ndvi, lon_min, lat_min, lon_max, lat_max)) in zip(axes, ndvi_data.items()):
#     img = ax.imshow(ndvi, extent=[lon_min, lon_max, lat_min, lat_max], cmap="RdYlGn", vmin=vmin, vmax=vmax)
#     ax.set_title(f"NDVI - {date}")
#     ax.set_xlabel("Longitude")
#     ax.set_ylabel("Latitude")

# # ✅ Final fix: Increase spacing on the right
# fig.subplots_adjust(right=0.8)  # Create more space on the right for the colorbar

# # ✅ Move the colorbar further right and reduce size
# cbar_ax = fig.add_axes([0.82, 0.3, 0.02, 0.4])  # Moves it right and shortens height
# fig.colorbar(img, cax=cbar_ax, label="NDVI Value")

# # ✅ SAVE THE FIGURE
# output_path = "/Users/hamza/Documents/NDVI_Comparison1.png"  # Change this to your preferred path
# plt.savefig(output_path, dpi=300, bbox_inches="tight")  # Save as high-quality image

# plt.show()

# print(f"✅ NDVI Comparison saved as: {output_path}")





# #NDVI Calculation 

# import rasterio
# import numpy as np

# # List of NDVI file paths
# ndvi_files = {
#     "January 2023": "/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jan2023.tif",
#     "July 2023": "/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jul2023.tif",
#     "January 2024": "/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jan2024.tif",
#     "July 2024": "/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jul2024.tif"
# }

# # Function to calculate NDVI statistics
# def calculate_ndvi_stats(file_path):
#     with rasterio.open(file_path) as src:
#         ndvi = src.read(1)  # Read NDVI band
#         ndvi = np.where((ndvi >= -1) & (ndvi <= 1), ndvi, np.nan)  # Remove invalid values

#     stats = {
#         "Min NDVI": np.nanmin(ndvi),
#         "Max NDVI": np.nanmax(ndvi),
#         "Mean NDVI": np.nanmean(ndvi),
#         "Median NDVI": np.nanmedian(ndvi),
#         "Std Dev": np.nanstd(ndvi)
#     }
#     return stats

# # Calculate and display results
# for date, path in ndvi_files.items():
#     stats = calculate_ndvi_stats(path)
#     print(f"\n📅 {date}")
#     for key, value in stats.items():
#         print(f"{key}: {value:.4f}")



# #NDVI_Difference_calculation

# import rasterio
# import numpy as np
# import matplotlib.pyplot as plt

# # Load two NDVI layers
# with rasterio.open("/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jan2024.tif") as src1:
#     ndvi_2024 = src1.read(1)
#     profile = src1.profile

# with rasterio.open("/Users/hamza/Documents/Thesis/NDVI/NDVI_Nile_Delta_Jan2023.tif") as src2:
#     ndvi_2023 = src2.read(1)

# # Calculate the NDVI difference
# ndvi_diff = ndvi_2024 - ndvi_2023

# # Plot the difference
# plt.figure(figsize=(10, 6))
# plt.imshow(ndvi_diff, cmap='bwr', vmin=-0.3, vmax=0.3)
# plt.colorbar(label="NDVI Change (2024 - 2023)")
# plt.title("NDVI Change: January 2024 - January 2023")
# plt.xlabel("Longitude")
# plt.ylabel("Latitude")
# plt.grid(False)
# plt.savefig("/Users/hamza/Documents/Thesis/NDVI/NDVI_Difference_Jan2024_Jan2023.png", dpi=300, bbox_inches='tight')
# plt.show()

# #Histogram
# plt.hist(ndvi_diff.flatten(), bins=50, color='Blue')
# plt.title("NDVI Change Distribution")
# plt.xlabel("NDVI Difference (2024 - 2023)")
# plt.ylabel("Frequency")
# plt.grid(False)
# plt.savefig("/Users/hamza/Documents/Thesis/NDVI/NDVI_Difference_Histo_Jan2024_Jan2023.png", dpi=300, bbox_inches='tight')
# plt.show()

# #Mean_NDVI-diff_value 
# mean_diff = np.nanmean(ndvi_diff)
# print(f"🟢 Mean NDVI Difference (2024 - 2023): {mean_diff:.4f}")



# #NDVI_5Years_chart

# # Load the Excel file
# import pandas as pd
# import matplotlib.pyplot as plt
# import datetime

# # Load the Excel file
# file_path = '/Users/hamza/Documents/Thesis/Comulative/MODIS_Monthly_NDVI_KafrElSheikh_2014_2024.xlsx'
# df = pd.read_excel(file_path, sheet_name='MODIS_Monthly_NDVI_KafrElSheikh')

# # 2. Rename and convert date column
# df = df.rename(columns={'Year': 'Date'})
# df['Date'] = pd.to_datetime(df['Date'])
# df['year'] = df['Date'].dt.year

# # 3. Align all years to the same year (2020) for x-axis comparison
# df['aligned_date'] = df['Date'].apply(lambda d: d.replace(year=2020))

# # 4. Plot setup
# plt.figure(figsize=(18, 8))
# line_styles = ['-', '--', '-.', ':']
# markers = ['o', 's', 'v', 'd', '^', 'P', 'X', '*', 'H', '>', '<', '.']
# colors = plt.cm.tab20.colors  # up to 20 different colors

# # 5. Plot each year's NDVI
# for idx, year in enumerate(sorted(df['year'].unique())):
#     df_year = df[df['year'] == year].sort_values('aligned_date')
#     plt.plot(
#         df_year['aligned_date'],
#         df_year['Mean_NDVI'],
#         label=str(year),
#         linestyle=line_styles[idx % len(line_styles)],
#         marker=markers[idx % len(markers)],
#         color=colors[idx % len(colors)],
#         linewidth=2,
#         markersize=5
#     )

# # 6. Customize x-axis to show all months twice (1st and 15th)
# tick_dates = []
# tick_labels = []
# for month in range(1, 13):
#     first = datetime.datetime(2020, month, 1)
#     mid = datetime.datetime(2020, month, 15)
#     tick_dates.extend([first, mid])
#     tick_labels.extend([first.strftime('%b %d'), mid.strftime('%b %d')])

# plt.xticks(tick_dates, tick_labels, rotation=45, fontsize=10)
# plt.yticks(fontsize=12)
# plt.xlabel("Month", fontsize=14)
# plt.ylabel("NDVI", fontsize=14)
# plt.title("Monthly NDVI Seasonal Comparison (2014–2024) – Kafr El-Sheikh", fontsize=16)
# plt.grid(True, linestyle='--', alpha=0.5)
# plt.legend(title='Year', fontsize=10, title_fontsize=11, loc='lower right', bbox_to_anchor=(0.95, 0.05))
# plt.tight_layout()
# # Save the figure before showing it
# plt.savefig(
#     '/Users/hamza/Documents/Thesis/Comulative/NDVI_Seasonal_Comparison_2014_2024.png',  # Change this path if needed
#     dpi=300,         # High quality for printing/pubs
#     bbox_inches='tight'  # Avoid cutting off labels
# )

# plt.show()





# #Future_Temperature
# import xarray as xr
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

# # === Load NetCDF Dataset ===
# ds = xr.open_dataset('/Users/hamza/Downloads/1893a1d7c209248f82461ecac67da3e2/tas_Amon_ACCESS-CM2_ssp585_r1i1p1f1_gn_20500116-20991216.nc')

# temp_k = ds['tas']  # Temperature in Kelvin
# temp_c = temp_k - 273.15  # Convert to Celsius

# # === Subset to Kafr El-Sheikh (nearest point) ===
# subset = temp_c.sel(lat=31.0625, lon=30.0625, method="nearest").squeeze()

# # === Time and temperature values ===
# time = pd.to_datetime(subset['time'].values)
# temps = subset.values

# # === Prepare plot ===
# fig, ax = plt.subplots(figsize=(14, 5))
# ax.plot(time, temps, color='darkred', linewidth=1.5, label="Monthly Temp")

# # === Annotate peak (summer) and trough (winter) per year ===
# df = pd.DataFrame({'time': time, 'temp': temps})
# df['year'] = df['time'].dt.year

# for year, group in df.groupby('year'):
#     max_row = group.loc[group['temp'].idxmax()]
#     min_row = group.loc[group['temp'].idxmin()]

#     # Annotate max (summer)
#     ax.annotate(f"{max_row['temp']:.1f}°C", 
#                 (max_row['time'], max_row['temp']),
#                 xytext=(0, 10), textcoords='offset points',
#                 ha='center', fontsize=5, color='darkgreen',
#                 bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='darkgreen'))

#     # Annotate min (winter)
#     ax.annotate(f"{min_row['temp']:.1f}°C", 
#                 (min_row['time'], min_row['temp']),
#                 xytext=(0, -15), textcoords='offset points',
#                 ha='center', fontsize=5, color='blue',
#                 bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='blue'))

# # === Labels and formatting ===
# ax.set_title("Monthly Avg Temp for Kafr El-Sheikh (2050–2099)", fontsize=14)
# ax.set_xlabel("Time", fontsize=12)
# ax.set_ylabel("Temperature (°C)", fontsize=12)
# ax.grid(True, linestyle='--', alpha=0.5)
# ax.legend(loc="upper left")

# # Improve x-axis ticks: show years only
# years = pd.date_range(start=time[0], end=time[-1], freq='YS')
# ax.set_xticks(years)
# ax.set_xticklabels([dt.strftime('%Y') for dt in years], rotation=45, fontsize=9)


# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/KafrElSheikh_MonthlyTemp_2050_2099.png", dpi=300)
# plt.show()



# #Future_annual_Trend
# import xarray as xr
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.stats import linregress

# # === Load NetCDF dataset ===
# ds = xr.open_dataset("/Users/hamza/Downloads/d4df10b331f84ab6bea5926bda7ee583/tas_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_19500116-19991216.nc")  # ← update path
# temp_k = ds['tas']  # Temperature in Kelvin
# temp_c = temp_k - 273.15  # Convert to Celsius

# # === Subset to Kafr El-Sheikh ===
# subset = temp_c.sel(lat=31.0625, lon=30.0625, method="nearest").squeeze()
# time = pd.to_datetime(subset['time'].values)
# temps = subset.values

# # === Create DataFrame and calculate yearly averages ===
# df = pd.DataFrame({'time': time, 'temp': temps})
# df['year'] = df['time'].dt.year
# yearly_avg = df.groupby('year')['temp'].mean().reset_index()

# # === Perform linear regression ===
# slope, intercept, r_value, p_value, std_err = linregress(yearly_avg['year'], yearly_avg['temp'])
# trend_line = intercept + slope * yearly_avg['year']

# # === Plot trend ===
# fig, ax = plt.subplots(figsize=(12, 5))
# ax.plot(yearly_avg['year'], yearly_avg['temp'], label='Yearly Avg Temp', color='darkred')
# ax.plot(yearly_avg['year'], trend_line, '--', color='navy', label=f"Trend: {slope*10:.2f} °C/decade")

# # === Labels & Aesthetics ===
# ax.set_title("Annual Avg Temp Trend – Kafr El-Sheikh (1950–1999)", fontsize=14)
# ax.set_xlabel("Year", fontsize=12)
# ax.set_ylabel("Temperature (°C)", fontsize=12)
# ax.grid(True, linestyle='--', alpha=0.5)
# ax.legend()
# plt.tight_layout()

# # === Save plot ===
# plt.savefig("/Users/hamza/Documents/Thesis/KafrElSheikh_TempTimeseries_1950_1999.png", dpi=300)
# plt.show()




# #Temperature_Timeseries_Copernicus
# import xarray as xr
# import pandas as pd
# import matplotlib.pyplot as plt

# # === Load Dataset ===
# ds = xr.open_dataset('//Users/hamza/Downloads/1893a1d7c209248f82461ecac67da3e2/tas_Amon_ACCESS-CM2_ssp585_r1i1p1f1_gn_20500116-20991216.nc')
# temp_k = ds['tas']
# temp_c = temp_k - 273.15

# # === Subset to Kafr El-Sheikh ===
# subset = temp_c.sel(lat=31.0625, lon=30.0625, method="nearest").squeeze()

# # === Time and temperature values ===
# time = pd.to_datetime(subset['time'].values)
# temps = subset.values

# # === Create dataframe and group by year ===
# df = pd.DataFrame({'time': time, 'temp': temps})
# df['year'] = df['time'].dt.year

# # === Extract yearly max and min ===
# summary = df.groupby('year').agg(Max_Temp=('temp', 'max'), Min_Temp=('temp', 'min')).reset_index()

# # === Export to CSV ===
# summary.to_csv("/Users/hamza/Documents/Thesis/KafrElSheikh_Yearly_Temp_Summary_2050_2099.csv", index=False)

# # === Plot Monthly Time Series (Clean Graph) ===
# fig, ax = plt.subplots(figsize=(14, 5))
# ax.plot(time, temps, color='darkred', linewidth=1.5, label="Monthly Temp")

# # Format
# ax.set_title("Monthly Avg Temp for Kafr El-Sheikh (1950-1999)", fontsize=14)
# ax.set_xlabel("Time", fontsize=12)
# ax.set_ylabel("Temperature (°C)", fontsize=12)
# ax.grid(True, linestyle='--', alpha=0.5)
# ax.legend(loc="upper left")

# # X-axis: yearly ticks
# years = pd.date_range(start=time[0], end=time[-1], freq='YS')
# ax.set_xticks(years)
# ax.set_xticklabels([dt.strftime('%Y') for dt in years], rotation=45, fontsize=9)

# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/KafrElSheikh_MonthlyTemp_2050.png", dpi=300)
# plt.show()




# #Monthly_Precipitation_Time_Series_for_Kafr_El-Sheikh
# import xarray as xr
# import pandas as pd
# import matplotlib.pyplot as plt

# # Load precipitation dataset
# ds = xr.open_dataset("/Users/hamza/Documents/Thesis/Copernicus/Prec-Historical/pr_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_19500116-19991216.nc")  # <- update this path
# precip_k = ds['pr']  # <- replace 'pr' with actual variable name

# # Convert from kg/m²/s to mm/month (approximate: 30 days/month)
# precip_mm_month = precip_k * 2.592e6  # 60 * 60 * 24 * 30

# # Subset to Kafr El-Sheikh
# subset = precip_mm_month.sel(lat=31.0625, lon=30.0625, method="nearest").squeeze()

# # Extract time and values
# time = pd.to_datetime(subset['time'].values)
# precips = subset.values

# # Plot monthly time series
# fig, ax = plt.subplots(figsize=(14, 5))
# ax.plot(time, precips, color='blue', linewidth=1.5, label="Monthly Precip")

# ax.set_title("Monthly Precipitation for Kafr El-Sheikh (1950–1999)", fontsize=14)
# ax.set_xlabel("Time", fontsize=12)
# ax.set_ylabel("Precipitation (mm/month)", fontsize=12)
# ax.grid(True, linestyle='--', alpha=0.5)
# ax.legend(loc='upper left')

# # Format x-axis ticks
# years = pd.date_range(start=time.min(), end=time.max(), freq='YS')
# ax.set_xticks(years)
# ax.set_xticklabels([dt.strftime('%Y') for dt in years], rotation=45, fontsize=9)

# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/KafrElSheikh_Monthly_Precip_1950-1999.png", dpi=300)
# plt.show()


# #Export_Prec_CSV_file
# import xarray as xr
# import pandas as pd

# # Load precipitation dataset
# ds = xr.open_dataset("/Users/hamza/Documents/Thesis/Copernicus/Prec-Historical/pr_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_19500116-19991216.nc")  # <- update this path
# precip_k = ds['pr']  # <- replace 'pr' with actual variable name

# # Convert from kg/m²/s to mm/month
# precip_mm_month = precip_k * 2.592e6

# # Subset to Kafr El-Sheikh
# subset = precip_mm_month.sel(lat=31.0625, lon=30.0625, method="nearest").squeeze()

# # Extract time and values
# time = pd.to_datetime(subset['time'].values)
# precips = subset.values

# # Create DataFrame
# df = pd.DataFrame({'time': time, 'precip': precips})
# df['year'] = df['time'].dt.year

# # Group by year and summarize
# summary = df.groupby('year').agg(
#     Total_Precip=('precip', 'sum'),
#     Max_Monthly=('precip', 'max'),
#     Min_Monthly=('precip', 'min')
# ).reset_index()

# # Export to CSV
# summary.to_csv("/Users/hamza/Documents/Thesis/KafrElSheikh_Yearly_Precip_Summary_1950_1999.csv", index=False)



# #Prec_Trend_line
# import xarray as xr
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.stats import linregress

# # Load NetCDF precipitation data
# ds = xr.open_dataset("/Users/hamza/Documents/Thesis/Copernicus/Prec-Historical/pr_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_19500116-19991216.nc")  # <- update path
# precip_k = ds['pr']  # <- update with actual variable name

# # Convert to mm/month (kg/m²/s × seconds per month)
# precip_mm_month = precip_k * 2.592e6

# # Subset to Kafr El-Sheikh
# subset = precip_mm_month.sel(lat=31.0625, lon=30.0625, method="nearest").squeeze()

# # Extract time and values
# time = pd.to_datetime(subset['time'].values)
# precips = subset.values

# # Create DataFrame
# df = pd.DataFrame({'time': time, 'precip': precips})
# df['year'] = df['time'].dt.year

# # Calculate annual total precipitation
# annual = df.groupby('year')['precip'].sum().reset_index()
# annual.columns = ['year', 'Total_Precip']

# # Perform linear regression
# slope, intercept, r_value, p_value, std_err = linregress(annual['year'], annual['Total_Precip'])
# annual['trend'] = intercept + slope * annual['year']

# # Plot
# fig, ax = plt.subplots(figsize=(10, 5))
# ax.plot(annual['year'], annual['Total_Precip'], marker='o', linestyle='-', label='Annual Total Precipitation', color='blue')
# ax.plot(annual['year'], annual['trend'], color='red', linestyle='--', linewidth=2, label=f'Trend (slope: {slope:.2f} mm/year)')

# # Format
# ax.set_title("Annual Precipitation Trend for Kafr El-Sheikh 1950-1999", fontsize=14)
# ax.set_xlabel("Year", fontsize=12)
# ax.set_ylabel("Total Precipitation (mm/year)", fontsize=12)
# ax.grid(True, linestyle='--', alpha=0.5)
# ax.legend()

# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/KafrElSheikh_Annual_Precip_Trend_1950-1999.png", dpi=300)
# plt.show()




# #ThreeModelsAnnualTemp
# import xarray as xr
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import linregress

# # === Coordinates for Kafr El-Sheikh ===
# target_lat = 31.1
# target_lon = 30.95

# # === Historical CMIP6 files ===
# model_files = {
#     "CNRM-CM6-1": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Temp/CNRM-CM6-1-2015/tas_Amon_CNRM-CM6-1_ssp126_r1i1p1f2_gr_20150116-20491216.nc",
#     "CanESM5": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Temp/CanESM5-2015/tas_Amon_CanESM5_ssp126_r1i1p1f1_gn_20150116-20491216.nc",
#     "CMCC-ESM2": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Temp/CMCC-ESM2-2015/tas_Amon_CMCC-ESM2_ssp126_r1i1p1f1_gn_20150116-20491216.nc"
# }

# model_dfs = []

# for model, file in model_files.items():
#     ds = xr.open_dataset(file, decode_times=False)

#     # Fix longitude if dataset is in 0–360 range
#     if ds.lon.max() > 180:
#         lon = target_lon + 360
#     else:
#         lon = target_lon

#     # Extract nearest grid point
#     ts = ds['tas'].sel(lat=target_lat, lon=lon, method='nearest')

#     # Assume monthly data: convert to yearly mean
#     n_months = len(ts.time)
#     n_years = n_months // 12
#     years = np.arange(2015, 2015 + n_years)
#     tas_vals = ts.values[:n_years * 12].reshape(n_years, 12)
#     annual_means = tas_vals.mean(axis=1)

#     df = pd.DataFrame({'year': years, model: annual_means})
#     model_dfs.append(df)

# # === Merge all models ===
# merged = model_dfs[0]
# for df in model_dfs[1:]:
#     merged = pd.merge(merged, df, on='year')

# # === Ensemble Mean and Trend ===
# merged['EnsembleMean'] = merged.iloc[:, 1:].mean(axis=1)
# slope, intercept, _, _, _ = linregress(merged['year'], merged['EnsembleMean'])
# merged['Trend'] = intercept + slope * merged['year']

# # === Convert from Kelvin to Celsius ===
# merged['EnsembleMean_C'] = merged['EnsembleMean'] - 273.15
# merged['Trend_C'] = merged['Trend'] - 273.15

# # === Save CSV ===
# merged.to_csv("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/TempAnnualtrend/Future/SSP2.6/2015-2049/KafrElSheikh_Tas_Ensemble_Historical_2015_2049-SSP2.6.csv", index=False)

# # === Plot Style Matching Your Example ===
# slope_decade = slope * 10  # °C per decade
# trend_label = f"Trend: {slope_decade:.2f} °C/decade"

# plt.figure(figsize=(12, 5))
# plt.plot(merged['year'], merged['EnsembleMean_C'], color='#8B0000', label='Yearly Avg Temp')
# plt.plot(merged['year'], merged['Trend_C'], 'b--', label=trend_label)
# plt.title("Annual Avg Temp Trend – Kafr El-Sheikh (2015–2049) SSP2.6 ", fontsize=14)
# plt.xlabel("Year", fontsize=12)
# plt.ylabel("Temperature (°C)", fontsize=12)
# plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.6)
# plt.legend(loc='upper left', fontsize=11)
# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/TempAnnualtrend/Future/SSP2.6/2015-2049/KafrElSheikh_Tas_Trend_Ensemble_2015_2049-SSP2.6.png", dpi=300)
# plt.show()



# #Monthly_3Models_Temp
# import xarray as xr
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

# # Coordinates for Kafr El-Sheikh
# target_lat = 31.1
# target_lon = 30.95

# # File paths for historical monthly tas data
# model_files = {
#     "CNRM-CM6-1": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Temp/CMCC-ESM2/tas_Amon_CMCC-ESM2_ssp126_r1i1p1f1_gn_20500116-20991216.nc",
#     "CanESM5": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Temp/CanESM5/tas_Amon_CanESM5_ssp126_r1i1p1f1_gn_20500116-20991216.nc",
#     "CMCC-ESM2": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Temp/CMCC-ESM2/tas_Amon_CMCC-ESM2_ssp126_r1i1p1f1_gn_20500116-20991216.nc"
# }

# monthly_df = pd.DataFrame()

# for model, path in model_files.items():
#     ds = xr.open_dataset(path, decode_times=False)

#     # Handle longitudes
#     lon_adj = target_lon + 360 if ds.lon.max() > 180 else target_lon

#     # Extract time series for nearest grid point
#     ts = ds['tas'].sel(lat=target_lat, lon=lon_adj, method='nearest')

#     # Extract monthly tas values and convert from K to °C
#     temps_C = ts.values - 273.15
#     monthly_df[model] = temps_C

# # Average across models
# monthly_df['EnsembleMean_C'] = monthly_df.mean(axis=1)

# # Create time index (monthly from Jan 1950 to Dec 1999)
# time_index = pd.date_range(start="2050-01", periods=len(monthly_df), freq='MS')
# monthly_df['time'] = time_index
# monthly_df.set_index('time', inplace=True)

# # Save to CSV if needed
# monthly_df.to_csv("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/MonthlyTemp/Future/SSP2.6/KafrElSheikh_MonthlyTemp_Ensemble_2050_2099.csv")

# # Plot
# plt.figure(figsize=(16, 5))
# plt.plot(monthly_df.index, monthly_df['EnsembleMean_C'], color='#8B0000', linewidth=1.5, label="Monthly Temp")
# plt.title("Monthly Avg Temp for Kafr El-Sheikh (2050–2099) SSP2.6", fontsize=14)
# plt.xlabel("Time", fontsize=12)
# plt.ylabel("Temperature (°C)", fontsize=12)
# plt.grid(True, linestyle='--', alpha=0.5)
# plt.legend(fontsize=11)
# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/MonthlyTemp/Future/SSP2.6/KafrElSheikh_MonthlyTemp_2050_2099.png", dpi=300)
# plt.show()






# #Hist_Prec_Annual_trend
# import xarray as xr
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import linregress

# # === Coordinates for Kafr El-Sheikh ===
# target_lat = 31.1
# target_lon = 30.95

# # === Model file paths (replace with your actual .nc paths) ===
# model_files = {
#     "CNRM-CM6-1": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Prec-2015/CNRM-CM6-1/pr_Amon_CNRM-CM6-1_ssp126_r1i1p1f2_gr_20150116-20491216.nc",
#     "CanESM5": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Prec-2015/CanESM5/pr_Amon_CanESM5_ssp126_r1i1p1f1_gn_20150116-20491216.nc",
#     "CMCC-ESM2": "/Users/hamza/Documents/Thesis/3Models_SSP2.6_Prec-2015/CMCC-ESM2/pr_Amon_CMCC-ESM2_ssp126_r1i1p1f1_gn_20150116-20491216.nc"
# }


# # === Step 1: Extract annual total precipitation in mm/year ===
# model_dfs = []

# for model, file in model_files.items():
#     ds = xr.open_dataset(file, decode_times=False)

#     # Adjust longitude if dataset is 0–360
#     lon = target_lon + 360 if ds.lon.max() > 180 else target_lon

#     # Extract nearest grid point
#     ts = ds['pr'].sel(lat=target_lat, lon=lon, method='nearest')

#     # Convert from kg/m²/s to mm/month
#     pr_vals = ts.values * 30 * 24 * 3600  # ≈ mm/month

#     # Reshape to (years, 12 months)
#     n_months = len(pr_vals)
#     n_years = n_months // 12
#     years = np.arange(2015, 2015 + n_years)
#     pr_vals = pr_vals[:n_years * 12].reshape(n_years, 12)

#     # Annual totals in mm
#     annual_pr = pr_vals.sum(axis=1)

#     # Store as DataFrame
#     df = pd.DataFrame({'year': years, model: annual_pr})
#     model_dfs.append(df)

# # === Step 2: Merge models & calculate ensemble ===
# merged = model_dfs[0]
# for df in model_dfs[1:]:
#     merged = pd.merge(merged, df, on='year')

# # Ensemble mean
# merged['EnsembleMean'] = merged.iloc[:, 1:].mean(axis=1)

# # Linear trend
# slope, intercept, _, _, _ = linregress(merged['year'], merged['EnsembleMean'])
# merged['Trend'] = intercept + slope * merged['year']

# # Save CSV
# merged.to_csv("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/PrecAnnualtrend/Future/SSP2.6/2015-2049/KafrElSheikh_Precip_Ensemble_SSP2.6_2015-2049.csv", index=False)

# # === Step 3: Styled Plot ===
# plt.figure(figsize=(12, 6))

# # Precipitation curve
# plt.plot(merged['year'], merged['EnsembleMean'], 'o-', color='blue', label='Annual Total Precipitation')

# # Trend line
# plt.plot(merged['year'], merged['Trend'], 'r--', linewidth=2, label=f"Trend (slope: {slope:.2f} mm/year)")

# # Titles and labels
# plt.title("Annual Precipitation Trend for Kafr El-Sheikh 2015-2049 - SSP2.6 ", fontsize=16)
# plt.xlabel("Year", fontsize=14)
# plt.ylabel("Total Precipitation (mm/year)", fontsize=14)

# # Grid and legend
# plt.grid(True, linestyle='--', alpha=0.7)
# plt.legend(fontsize=12)
# plt.tight_layout()

# # Save + show
# plt.savefig("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/PrecAnnualtrend/Future/SSP2.6/2015-2049/KafrElSheikh_Precip_Trend_SSP2.6_2015_2049.png", dpi=300)
# plt.show()




# #Monthly_3Models_Prec
# import xarray as xr
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

# # === Coordinates for Kafr El-Sheikh ===
# target_lat = 31.1
# target_lon = 30.95

# # === Paths to model files ===
# model_files = {
#     "CNRM-CM6-1": "/Users/hamza/Documents/Thesis/3Models_SSP8.5_Prec-2049/CNRM-CM6-1 /pr_Amon_CNRM-CM6-1_ssp585_r1i1p1f2_gr_20150116-20491216.nc",
#     "CanESM5": "/Users/hamza/Documents/Thesis/3Models_SSP8.5_Prec-2049/CanESM5/pr_Amon_CanESM5_ssp585_r1i1p1f1_gn_20150116-20491216.nc",
#     "CMCC-ESM2": "/Users/hamza/Documents/Thesis/3Models_SSP8.5_Prec-2049/CMCC-ESM2/pr_Amon_CMCC-ESM2_ssp585_r1i1p1f1_gn_20150116-20491216.nc"
# }

# # === Extract monthly precipitation (in mm/month) ===
# monthly_df = pd.DataFrame()

# for model, path in model_files.items():
#     ds = xr.open_dataset(path, decode_times=False)

#     # Fix longitude
#     lon = target_lon + 360 if ds.lon.max() > 180 else target_lon

#     # Extract nearest grid cell
#     ts = ds['pr'].sel(lat=target_lat, lon=lon, method='nearest')

#     # Convert to mm/month (kg/m²/s × seconds in ~30 days)
#     pr_mm = ts.values * 30 * 24 * 3600
#     monthly_df[model] = pr_mm

# # === Average the models ===
# monthly_df['EnsembleMean'] = monthly_df.mean(axis=1)

# # === Create monthly datetime index ===
# time_index = pd.date_range(start='2015-01', periods=len(monthly_df), freq='MS')
# monthly_df['Time'] = time_index
# monthly_df.set_index('Time', inplace=True)

# # === Save to CSV ===
# monthly_df.to_csv("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/MonthlyPrec/Future/SSP8.5/KafrElSheikh_Precip_Monthly_SSP8.5_2015_2049.csv")

# # === Plot ===
# plt.figure(figsize=(18, 5))
# plt.plot(monthly_df.index, monthly_df['EnsembleMean'], color='blue', linewidth=1, label='Monthly Precip')
# plt.title("Monthly Precipitation for Kafr El-Sheikh (2015–2049) SSP8.5", fontsize=15)
# plt.xlabel("Time", fontsize=12)
# plt.ylabel("Precipitation (mm/month)", fontsize=12)
# plt.legend(fontsize=11)
# plt.grid(True, linestyle='--', alpha=0.7)
# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/MonthlyPrec/Future/SSP8.5/KafrElSheikh_Monthly_Precip_SSP8.5_2015-2049.png", dpi=300)
# plt.show()



# #Avg_Temp_Annual_comparison_graph
# import pandas as pd
# import matplotlib.pyplot as plt

# # === Load formatted CSVs ===
# files = {
#     # "1950–1999 Historical": pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_Historical_1950_1999.csv"),
#     # "2000–2014 Historical":  pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_Historical_2000_2014.csv"),
#     "2000–2049 SSP1-2.6":  pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_SSP126_1950_2049.csv"),
#     "2000–2049 SSP5-8.5":  pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_SSP585_1950_2049.csv"),
#     "2050–2099 SSP1-2.6":  pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_SSP126_2050_2099.csv"),
#     "2050–2099 SSP5-8.5":  pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_SSP585_2050_2099.csv")
# }

# # === Combine all data ===
# all_data = pd.DataFrame()
# for label, df in files.items():
#     df['label'] = label
#     all_data = pd.concat([all_data, df], ignore_index=True)


# # === Plotting ===
# plt.figure(figsize=(14, 6))

# # Custom color map
# colors = {
#     "1950–1999 Historical": "#444444",
#     "2000–2049 SSP1-2.6": "#3d9257",
#     "2000–2049 SSP5-8.5": "#aa1e2d",
#     "2050–2099 SSP1-2.6": "#65c07a",
#     "2050–2099 SSP5-8.5": "#d65454"
# }

# # Plot each label group
# for label in all_data['label'].unique():
#     subset = all_data[all_data['label'] == label]
#     plt.plot(subset['year'], subset['EnsembleMean_C'], label=label, color=colors[label], linewidth=2)

# # === Final plot formatting ===
# plt.title("Annual Avg Temperature in Kafr El-Sheikh (Historical & CMIP6 SSPs)", fontsize=15)
# plt.xlabel("Year", fontsize=12)
# plt.ylabel("Temperature (°C)", fontsize=12)
# plt.grid(True, linestyle='--', alpha=0.6)
# plt.legend(title="Scenario", fontsize=10)
# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/KafrElSheikh_Temp_Comparison_Historical_SSPs.png", dpi=300)
# plt.show()






# #Avg_Temp_Annual_comparison_graph
# import pandas as pd
# import matplotlib.pyplot as plt

# # Load & combine by scenario
# historical = pd.concat([
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_Historical_1950_1999.csv"),
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_Historical_2000_2014.csv")
# ], ignore_index=True)
# historical['Scenario'] = 'Historical'


# ssp126 = pd.concat([
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_SSP126_2000_2049.csv"),
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_SSP126_2050_2099.csv")
# ], ignore_index=True)

# # 🛑 Only keep years 2015–2099
# ssp126["year"] = ssp126["year"].astype(int)
# ssp126 = ssp126[ssp126["year"] >= 2014]
# ssp126['Scenario'] = 'SSP1-2.6'



# ssp585 = pd.concat([
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_SSP585_2000_2049.csv"),
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_SSP585_2050_2099.csv")
# ], ignore_index=True)

# ssp585["year"] = ssp585["year"].astype(int)
# ssp585 = ssp585[ssp585["year"] >= 2014]
# ssp585['Scenario'] = 'SSP5-8.5'


# # Combine all
# df_all = pd.concat([historical, ssp126, ssp585], ignore_index=True)


# # Plot
# plt.figure(figsize=(14,6))
# colors = {
#     'Historical': '#333333',
#     'SSP1-2.6': '#2ca02c',
#     'SSP5-8.5': '#d62728'
# }

# for scenario, group in df_all.groupby('Scenario'):
#     plt.plot(group['year'], group['EnsembleMean_C'], label=scenario, color=colors[scenario], linewidth=2)

# plt.xlabel("Year", fontsize=12)
# plt.ylabel("Temperature (°C)", fontsize=12)
# plt.title("CMIP6 Annual Avg Temperature – Kafr El-Sheikh (1950–2099)", fontsize=14)
# plt.legend(title="Scenario", fontsize=11)
# plt.grid(True, linestyle='--', alpha=0.5)
# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/KafrElSheikh_Temp_Comparison_Historical_SSPs.png", dpi=300)
# plt.show()





# # Annual Precipitation Comparison Graph
# import pandas as pd
# import matplotlib.pyplot as plt

# # === Load cleaned files ===
# files = {
#     "Historical": [
#         pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_Formatted_Historical_1950_1999.csv"),
#         pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_Formatted_Historical_2000_2014.csv")
#     ],
#     "SSP1-2.6": [
#         pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_SSP126_2000_2049.csv").query("year >= 2014"),
#         pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_Formatted_SSP126_2050_2099.csv")
#     ],
#     "SSP5-8.5": [
#         pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_SSP585_2000_2049.csv").query("year >= 2014"),
#         pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_Formatted_SSP585_2050_2099.csv")
#     ]
# }



# ssp126 = pd.concat([
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_SSP126_2000_2049.csv"),
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_Formatted_SSP126_2050_2099.csv")
# ], ignore_index=True)

# ssp126["year"] = ssp126["year"].astype(int)
# ssp126 = ssp126[ssp126["year"] >= 2015]  # ✅ Only keep 2015–2099
# ssp126["Scenario"] = "SSP1-2.6"


# ssp585 = pd.concat([
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_SSP585_2000_2049.csv"),
#     pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Precip_Formatted_SSP585_2050_2099.csv")
# ], ignore_index=True)

# ssp585["year"] = ssp585["year"].astype(int)
# ssp585 = ssp585[ssp585["year"] >= 2015]  # ✅ Only keep 2015–2099
# ssp585["Scenario"] = "SSP5-8.5"

# # === Combine and tag ===
# all_data = pd.DataFrame()
# for scenario, dfs in files.items():
#     combined = pd.concat(dfs)
#     combined["Scenario"] = scenario
#     all_data = pd.concat([all_data, combined])

# # === Plot: clean scenario lines ===
# plt.figure(figsize=(14, 6))

# colors = {
#     "Historical": "black",
#     "SSP1-2.6": "green",
#     "SSP5-8.5": "red"
# }

# for scenario in all_data["Scenario"].unique():
#     df = all_data[all_data["Scenario"] == scenario]
#     df["EnsembleMean_mm"] = pd.to_numeric(df["EnsembleMean_mm"], errors="coerce")
#     df = df.groupby("year", as_index=False)["EnsembleMean_mm"].mean()
#     plt.plot(df["year"], df["EnsembleMean_mm"], label=scenario, color=colors[scenario], linewidth=2)

# # === Final formatting ===
# plt.title("Annual Precipitation in Kafr El-Sheikh (Historical & CMIP6 SSPs (1950–2099)", fontsize=15)
# plt.xlabel("Year", fontsize=12)
# plt.ylabel("Total Precipitation (mm/year)", fontsize=12)
# plt.grid(True, linestyle="--", alpha=0.5)
# plt.legend(title="Scenario", fontsize=10)
# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/KafrElSheikh_Precip_Comparison_Historical_SSPs.png", dpi=300)
# plt.show()


# #Mann-Kendall_Trend
# import pandas as pd
# import pymannkendall as mk
# import matplotlib.pyplot as plt

# # === 1. Load your temperature data ===
# df = pd.read_csv("/Users/hamza/Downloads/KafrElSheikh_Tas_Ensemble_Historical_1950_1999.csv")  # 🔁 Replace with your actual file path
# df = df.rename(columns={'year': 'Year', 'EnsembleMean_C': 'Temperature'})
# df = df.sort_values(by='Year')

# # === 2. Run the Mann-Kendall Test ===
# result = mk.original_test(df['Temperature'])

# # Extract results
# trend = result.trend
# tau = result.Tau
# p_value = result.p
# slope = result.slope
# h = result.h

# # === 3. Plot the time series ===
# plt.figure(figsize=(10, 5))
# plt.plot(df['Year'], df['Temperature'], marker='o', linestyle='-', label='Temperature')

# # === 4. Plot Sen’s Slope Trend Line ===
# start_temp = df['Temperature'].iloc[0]
# sen_line = start_temp + slope * (df['Year'] - df['Year'].iloc[0])
# plt.plot(df['Year'], sen_line, 'r--', label="Sen's Slope Trend")

# # === 5. Annotate Mann-Kendall Results ===
# textstr = '\n'.join((
#     f'Trend: {trend}',
#     f'Tau = {tau:.3f}',
#     f'p-value = {p_value:.3e}',
#     f"Slope = {slope:.3f} °C/year",
#     f'Significant: {"Yes" if h == 1 else "No"}'
# ))
# plt.text(0.02, 0.95, textstr, transform=plt.gca().transAxes,
#          fontsize=8, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))

# # === 6. Final Plot Styling ===
# plt.title('Mann-Kendall Trend in Annual Temperature - Hist (1950 - 1999)')
# plt.xlabel('Year')
# plt.ylabel('Temperature (°C)')
# plt.legend()
# plt.tight_layout()

# # === 7. Save the Plot ===
# plt.savefig("/Users/hamza/Documents/Thesis/Copernicus/3Models_results/Mann-Kendall results/Mann-Kendalltemperature_trend Hist (1950 - 1999).png", dpi=300)
# plt.show()

# print("✅ Plot saved as 'temperature_trend.png'")


# # RiceYield_Correlation
# import pandas as pd
# from scipy.stats import pearsonr, spearmanr
# import matplotlib.pyplot as plt

# # === Load your data ===
# df = pd.read_csv("/Users/hamza/Downloads/SSP585_Climate_Yield_Cleaned.csv")

# # === Calculate Correlation Coefficients ===
# pearson_temp = pearsonr(df["Kafr_Yield_kg_ha_Estimated"], df["Temperature"])
# spearman_temp = spearmanr(df["Kafr_Yield_kg_ha_Estimated"], df["Temperature"])

# pearson_precip = pearsonr(df["Kafr_Yield_kg_ha_Estimated"], df["Precipitation"])
# spearman_precip = spearmanr(df["Kafr_Yield_kg_ha_Estimated"], df["Precipitation"])

# # === Print Results ===
# print("✅ Pearson  (Yield vs Temp):", round(pearson_temp[0], 3))
# print("✅ Spearman (Yield vs Temp):", round(spearman_temp.correlation, 3))
# print("✅ Pearson  (Yield vs Precip):", round(pearson_precip[0], 3))
# print("✅ Spearman (Yield vs Precip):", round(spearman_precip.correlation, 3))

# # === Optional: Plot the relationships ===
# plt.figure(figsize=(14, 5))

# # Yield vs Temperature
# plt.subplot(1, 2, 1)
# plt.scatter(df["Temperature"], df["Kafr_Yield_kg_ha_Estimated"], c="red")
# plt.title("Yield vs Temperature")
# plt.xlabel("Temperature (°C)")
# plt.ylabel("Yield (kg/ha)")

# # Yield vs Precipitation
# plt.subplot(1, 2, 2)
# plt.scatter(df["Precipitation"], df["Kafr_Yield_kg_ha_Estimated"], c="blue")
# plt.title("Yield vs Precipitation")
# plt.xlabel("Precipitation (mm)")
# plt.ylabel("Yield (kg/ha)")

# plt.tight_layout()
# plt.show()


# RiceYield_Correlation_Temperature
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import numpy as np

# Load the cleaned dataset
df = pd.read_csv("/Users/hamza/Documents/Thesis/Pearson and Spearman/SSP585_Rice_Yield.csv")

# Define column names
yield_col = "Kafr_Yield_kg_ha_Estimated"
temp_col = "Temperature"
precip_col = "Precipitation"

# Drop missing rows
df.dropna(subset=[yield_col, temp_col, precip_col], inplace=True)

# === Correlation Coefficients ===
pearson_temp_r, _ = pearsonr(df[temp_col], df[yield_col])
spearman_temp_r, _ = spearmanr(df[temp_col], df[yield_col])

pearson_precip_r, _ = pearsonr(df[precip_col], df[yield_col])
spearman_precip_r, _ = spearmanr(df[precip_col], df[yield_col])

# === Plotting ===
plt.figure(figsize=(14, 6))

# 1. Pearson (Linear): Temperature vs Yield
plt.subplot(1, 2, 1)
sns.regplot(x=temp_col, y=yield_col, data=df, color="orange", line_kws={"color": "black"})
plt.title(f"SSP5–8.5 | Pearson (Linear): r = {round(pearson_temp_r, 3)}")
plt.xlabel("Temperature (°C)")
plt.ylabel("Rice Yield (kg/ha)")

# 2. Spearman (Non-linear): Temperature vs Yield
plt.subplot(1, 2, 2)
x = df[temp_col].values.reshape(-1, 1)
y = df[yield_col].values.reshape(-1, 1)

poly = PolynomialFeatures(degree=2)
x_poly = poly.fit_transform(x)

model = LinearRegression().fit(x_poly, y)
x_range = np.linspace(x.min(), x.max(), 300).reshape(-1, 1)
y_pred = model.predict(poly.transform(x_range))

plt.scatter(x, y, color="red", alpha=0.6)
plt.plot(x_range, y_pred, color="black", linewidth=2)
plt.title(f"SSP5–8.5 | Spearman (Non-linear): ρ = {round(spearman_temp_r, 3)}")
plt.xlabel("Temperature (°C)")
plt.ylabel("Rice Yield (kg/ha)")

plt.tight_layout()
plt.savefig("/Users/hamza/Documents/Thesis/Pearson and Spearman/Temperature_Yield_Correlation_Plots-SSP8.5.png")
plt.show()




# # RiceYield_Correlation_Precipitation
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import pearsonr, spearmanr
# from sklearn.linear_model import LinearRegression
# from sklearn.preprocessing import PolynomialFeatures
# import numpy as np

# # Load the cleaned dataset
# df = pd.read_csv("/Users/hamza/Documents/Thesis/Pearson and Spearman/SSP585_Rice_Yield.csv")

# # Define column names
# yield_col = "Kafr_Yield_kg_ha_Estimated"
# precip_col = "Precipitation"

# # Drop missing values
# df.dropna(subset=[yield_col, precip_col], inplace=True)

# # Correlation Coefficients
# pearson_precip_r = pearsonr(df[precip_col], df[yield_col])
# spearman_precip_r = spearmanr(df[precip_col], df[yield_col])

# # Plotting
# plt.figure(figsize=(12, 5))

# # 1. Pearson: Precipitation vs Yield (Linear)
# plt.subplot(1, 2, 1)
# sns.regplot(x=precip_col, y=yield_col, data=df, color="blue", line_kws={"color": "black"})
# plt.title(f"SSP5–8.5 | Pearson (Linear): r = {round(pearson_precip_r[0], 3)}")
# plt.xlabel("Precipitation (mm)")
# plt.ylabel("Rice Yield (kg/ha)")

# # 2. Spearman: Precipitation vs Yield (Non-linear)
# x = df[precip_col].values.reshape(-1, 1)
# y = df[yield_col].values.reshape(-1, 1)

# poly = PolynomialFeatures(degree=2)
# x_poly = poly.fit_transform(x)

# model = LinearRegression().fit(x_poly, y)
# x_range = np.linspace(x.min(), x.max(), 300).reshape(-1, 1)
# x_range_poly = poly.transform(x_range)
# y_pred = model.predict(x_range_poly)

# plt.subplot(1, 2, 2)
# plt.scatter(x, y, color="deepskyblue", alpha=0.7)
# plt.plot(x_range, y_pred, color="black", linewidth=2)
# plt.title(f"SSP5–8.5 | Spearman (Non-linear): ρ = {round(spearman_precip_r.correlation, 3)}")
# plt.xlabel("Precipitation (mm)")
# plt.ylabel("Rice Yield (kg/ha)")

# plt.tight_layout()
# plt.savefig("/Users/hamza/Documents/Thesis/Pearson and Spearman/Precipitation_Yield_Correlation_Plots_SSP8.5.png")
# plt.show()









