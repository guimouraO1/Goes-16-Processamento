#------------------------------------------------------------------------------------------------------
# Required modules
import matplotlib.pyplot as plt              # Import the Matplotlib package
import numpy as np                           # Import the Numpy package
from netCDF4 import Dataset                  # Import the NetCDF Python interface
import cartopy, cartopy.crs as ccrs          # Plot maps
import cartopy.io.shapereader as shpreader   # Import shapefiles
import math                                  # Import math
from datetime import datetime, timedelta     # Library to convert julian day to dd-mm-yyyy
from truecolor import area_para_recorte
from remap import remap
import time as time    



print('Script started.')
start = time.time()  

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# File to read
path_ch08 = "OR_ABI-L2-CMIPF-M6C08_G16_s20232641020207_e20232641029515_c20232641029592.nc"
path_ch10 = "OR_ABI-L2-CMIPF-M6C10_G16_s20232641020207_e20232641029526_c20232641029583.nc"
path_ch12 = "OR_ABI-L2-CMIPF-M6C12_G16_s20232641020207_e20232641029521_c20232641029598.nc"
path_ch13 = "OR_ABI-L2-CMIPF-M6C13_G16_s20232641020207_e20232641029526_c20232641029588.nc"


# Read the file using the NetCDF library
file_ch08 = Dataset(path_ch08)
file_ch10 = Dataset(path_ch10)
file_ch12 = Dataset(path_ch12)
file_ch13 = Dataset(path_ch13)

variable = "CMI"
v_extent = 'br'
extent, resolution = area_para_recorte(v_extent)

# Read the resolution
band_resolution_km = getattr(file_ch08, 'spatial_resolution')
band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

# Division factor to reduce image size
f = math.ceil(float(resolution / band_resolution_km))

# Read the central longitude
longitude = file_ch08.variables['goes_imager_projection'].longitude_of_projection_origin
# Read the semi major axis
a = file_ch08.variables['goes_imager_projection'].semi_major_axis
# Read the semi minor axis
b = file_ch08.variables['goes_imager_projection'].semi_minor_axis

# Getting the file time and date
add_seconds = int(file_ch08.variables['time_bounds'][0])
date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
date_formated = date.strftime('%Y-%m-%d %H:%M UTC')
date_file = date.strftime('%Y%m%d%H%M')
year = date.strftime('%Y')
month = date.strftime('%m')
day = date.strftime('%d')
hour = date.strftime('%H')
minutes = date.strftime('%M')

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Call the reprojection funcion
grid = remap(path_ch08, variable, extent, resolution)
# Read the data returned by the function 
data_ch08 = grid.ReadAsArray() - 273.15
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Call the reprojection funcion
grid = remap(path_ch10, variable, extent, resolution)
# Read the data returned by the function 
data_ch10 = grid.ReadAsArray() - 273.15
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Call the reprojection funcion
grid = remap(path_ch12, variable, extent, resolution)
# Read the data returned by the function 
data_ch12 = grid.ReadAsArray() - 273.15
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Call the reprojection funcion
grid = remap(path_ch13, variable, extent, resolution)
# Read the data returned by the function 
data_ch13 = grid.ReadAsArray() - 273.15
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------


# RGB Components
R = data_ch08 - data_ch10
G = data_ch12 - data_ch13
B = data_ch08

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
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Plot configuration
plot_config = {
"resolution": band_resolution_km, 
"dpi": 150,  
"title_text": "GOES-16 AIRMASS RGB", "title_size": int(data_ch08.shape[1] * 0.005), "title_x_offset": int(data_ch08.shape[1] * 0.01), "title_y_offset": data_ch08.shape[0] - int(data_ch08.shape[0] * 0.016), 
"file_name_id_1": "G16",  "file_name_id_2": "ARMRGB" 
}
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

d_p_i = 150
fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

# Define the projection
proj = ccrs.PlateCarree()

# Use the PlateCarree projection in cartopy
ax = plt.axes([0, 0, 1, 1], projection=proj)

# Define the image extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]

# Plot the image
img = ax.imshow(RGB, origin='upper', extent=img_extent, zorder=3)
  
# Add a title
plt.annotate(plot_config["title_text"] + " " + date_formated , xy=(plot_config["title_x_offset"], plot_config["title_y_offset"]), xycoords='figure pixels', fontsize=plot_config["title_size"], fontweight='bold', color='white', bbox=dict(boxstyle="round",fc=(0.0, 0.0, 0.0), ec=(1., 1., 1.)), zorder=7)

# Save the image
plt.savefig(plot_config["file_name_id_1"] + "_" + plot_config["file_name_id_2"] + "_" + date_file + '.png', bbox_inches='tight', pad_inches=0, facecolor='black')