#######################################################################################################
# LICENSE
# Copyright (C) 2019 - INPE - NATIONAL INSTITUTE FOR SPACE RESEARCH - BRAZIL
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU 
# General Public License as published by the Free Software Foundation, either version 3 of the License, 
# or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
# Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. 
# If not, see http://www.gnu.org/licenses/.
#######################################################################################################
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Python Script: Airmass RGB
# Quick Guide: http://rammb.cira.colostate.edu/training/visit/quick_guides/QuickGuide_GOESR_AirMassRGB_final.pdf
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Required modules
import matplotlib.pyplot as plt              # Import the Matplotlib package
import numpy as np                           # Import the Numpy package
from netCDF4 import Dataset                  # Import the NetCDF Python interface
import cartopy, cartopy.crs as ccrs          # Plot maps
import cartopy.io.shapereader as shpreader   # Import shapefiles
import math                                  # Import math
from datetime import datetime, timedelta     # Library to convert julian day to dd-mm-yyyy
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# File to read
image1 = "OR_ABI-L2-CMIPF-M6C08_G16_s20191841700283_e20191841709591_c20191841710070.nc"

# Read the file using the NetCDF library
file1 = Dataset(image1)

# Desired resolution
resolution = 8

# Read the resolution
band_resolution_km = getattr(file1, 'spatial_resolution')
band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

# Division factor to reduce image size
f = math.ceil(float(resolution / band_resolution_km))

# Read the central longitude
longitude = file1.variables['goes_imager_projection'].longitude_of_projection_origin

# Calculate the image extent 
h = file1.variables['goes_imager_projection'].perspective_point_height
x = file1.variables['x_image_bounds'] * h 
y = file1.variables['y_image_bounds'] * h 

# Reading the file time and date
add_seconds = int(file1.variables['time_bounds'][0])
date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
date_formated = date.strftime('%Y-%m-%d %H:%M UTC')
date_file = date.strftime('%Y%m%d%H%M')
year = date.strftime('%Y')
month = date.strftime('%m')
day = date.strftime('%d')
hour = date.strftime('%H')
minutes = date.strftime('%M')

# Convert to Celsius
data1 = file1.variables['CMI'][:,:][::f ,::f] - 273.15
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# File to read
image2 = "OR_ABI-L2-CMIPF-M6C10_G16_s20191841700283_e20191841710002_c20191841710078.nc"

# Read the file using the NetCDF library
file2 = Dataset(image2)

# Convert to Celsius
data2 = file2.variables['CMI'][:,:][::f ,::f] - 273.15
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# File to read
image3 = "OR_ABI-L2-CMIPF-M6C12_G16_s20191841700283_e20191841709596_c20191841710078.nc"

# Read the file using the NetCDF library
file3 = Dataset(image3)

# Convert to Celsius
data3 = file3.variables['CMI'][:,:][::f ,::f] - 273.15
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# File to read
image4 = "OR_ABI-L2-CMIPF-M6C13_G16_s20191841700283_e20191841710002_c20191841710081.nc"

# Read the file using the NetCDF library
file4 = Dataset(image4)

# Convert to Celsius
data4 = file4.variables['CMI'][:,:][::f ,::f] - 273.15
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# RGB Components
R = data1 - data2
G = data3 - data4
B = data1

# Minimuns and Maximuns
Rmin = -26.2
Rmax = 0.6

Gmin = -43.2
Gmax = 6.7

Bmin = -29.25
Bmax = -64.65

R[R<Rmin] = Rmin
R[R>Rmax] = Rmax

G[G<Gmin] = Gmin
G[G>Gmax] = Gmax

B[B<Bmax] = Bmax
B[B>Bmin] = Bmin

# Choose the gamma
gamma = 1

# Normalize the data
R = ((R - Rmin) / (Rmax - Rmin)) ** (1/gamma)
G = ((G - Gmin) / (Gmax - Gmin)) ** (1/gamma)
B = ((B - Bmin) / (Bmax - Bmin)) ** (1/gamma) 

# Create the RGB
RGB = np.stack([R, G, B], axis=2)

# Eliminate values outside the globe
mask = (RGB == [R[0,0],G[0,0],B[0,0]]).all(axis=2)
RGB[mask] = np.nan
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Plot configuration
plot_config = {
"resolution": band_resolution_km, 
"dpi": 150, 
"states_color": 'grey', "states_width": data1.shape[0] * 0.00006, 
"countries_color": 'white', "countries_width": data1.shape[0] * 0.00012,
"continents_color": 'white', "continents_width": data1.shape[0] * 0.00025,
"grid_color": 'white', "grid_width": data1.shape[0] * 0.00025, "grid_interval": 10.0,
"title_text": "GOES-16 AIRMASS RGB", "title_size": int(data1.shape[1] * 0.005), "title_x_offset": int(data1.shape[1] * 0.01), "title_y_offset": data1.shape[0] - int(data1.shape[0] * 0.016), 
"file_name_id_1": "G16",  "file_name_id_2": "ARMRGB" 
}
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
fig = plt.figure(figsize=(data1.shape[1]/float(plot_config["dpi"]), data1.shape[0]/float(plot_config["dpi"])), dpi=plot_config["dpi"])

# Define the projection
proj = ccrs.Geostationary(central_longitude=longitude, satellite_height=h)
img_extent = (x.min(), x.max(), y.min(), y.max())

# Use the Geostationary projection in cartopy
ax = plt.axes([0, 0, 1, 1], projection=proj)

# Plot the image
img = ax.imshow(RGB, origin='upper', extent=img_extent, zorder=3)

# Add countries
shapefile = list(shpreader.Reader('ne_50m_admin_0_countries.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=plot_config["countries_color"],facecolor='none', linewidth=plot_config["countries_width"], zorder=4)

# Add continents
shapefile = list(shpreader.Reader('ne_10m_coastline.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=plot_config["continents_color"],facecolor='none', linewidth=plot_config["continents_width"], zorder=5)
  
# Add coastlines, borders and gridlines
ax.gridlines(color=plot_config["grid_color"], alpha=0.5, linestyle='--', linewidth=plot_config["grid_width"], xlocs=np.arange(-180, 180, plot_config["grid_interval"]), ylocs=np.arange(-180, 180, plot_config["grid_interval"]), draw_labels=False, zorder=6)

# Remove the outline border
ax.outline_patch.set_visible(False)
  
# Add a title
plt.annotate(plot_config["title_text"] + " " + date_formated , xy=(plot_config["title_x_offset"], plot_config["title_y_offset"]), xycoords='figure pixels', fontsize=plot_config["title_size"], fontweight='bold', color='white', bbox=dict(boxstyle="round",fc=(0.0, 0.0, 0.0), ec=(1., 1., 1.)), zorder=7)

# Save the image
plt.savefig(plot_config["file_name_id_1"] + "_" + plot_config["file_name_id_2"] + "_" + date_file + '.png', bbox_inches='tight', pad_inches=0, facecolor='black')
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------