import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset                                  # Read / Write NetCDF4 files
from mpl_toolkits.axes_grid1.inset_locator import inset_axes # Add a child inset axes to this existing axes.
from datetime import datetime, timedelta                     # Library to convert julian day to dd-mm-yyyy
from utilities import loadCPT                              # Import the CPT convert function
from matplotlib.colors import LinearSegmentedColormap        # Linear interpolation for color maps
import matplotlib.pyplot as plt                              # Plotting library
import matplotlib.colors                                     # Matplotlib colors
import numpy as np                                           # Scientific computing with Python
import cartopy, cartopy.crs as ccrs                          # Plot maps
import cartopy.io.shapereader as shpreader                   # Import shapefiles
import time as t                                             # Time access and conversion
import math                                                  # Import math
import re                                                    # re
import glob                                                  # Unix style pathname pattern expansion
import sys                                                   # Import the "system specific parameters and functions" module
import os 												     # Miscellaneous operating system interfaces
from os.path import dirname, abspath                         # Return a normalized absolutized version of the pathname path 
import platform                                              # To check which OS is being used                          # Update the HTML animation 
from remap import remap                                      # Import the Remap function
import warnings
warnings.filterwarnings("ignore")


print('Script started.')
start = t.time()  

dir_in = f'/mnt/e/truecolor/'
dir_shapefiles = f'{dir_in}shapefiles/'
dir_colortables = f'{dir_in}colortables/'
dir_logos = f'{dir_in}logos/'
dir_out = f'{dir_in}output/'

# Bands files
path_ch01 = f'{dir_in}CG_ABI-L2-CMIPF-M6C01_G16_s20232631130206_e20232631139514_c20232631139581.nc'
path_ch02 = f'{dir_in}CG_ABI-L2-CMIPF-M6C02_G16_s20232631130206_e20232631139514_c20232631139575.nc'
path_ch03 = f'{dir_in}CG_ABI-L2-CMIPF-M6C03_G16_s20232631130206_e20232631139514_c20232631139581.nc'

# For the log
path = path_ch01

print(path_ch01)
print(path_ch02)
print(path_ch03)

# Read the image
file_ch01 = Dataset(path_ch01)

# Read the satellite 
satellite = getattr(file_ch01, 'platform_ID')

# Read the band number
band = str(file_ch01.variables['band_id'][0]).zfill(2)

# Desired resolution
resolution = 3.0

# Read the resolution
band_resolution_km = getattr(file_ch01, 'spatial_resolution')
band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

# Division factor to reduce image size
f = math.ceil(float(resolution / band_resolution_km))

# Read the central longitude
longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin

# Read the semi major axis
a = file_ch01.variables['goes_imager_projection'].semi_major_axis

# Read the semi minor axis
b = file_ch01.variables['goes_imager_projection'].semi_minor_axis

# Calculate the image extent 
h = file_ch01.variables['goes_imager_projection'].perspective_point_height
x1 = file_ch01.variables['x_image_bounds'][0] * h 
x2 = file_ch01.variables['x_image_bounds'][1] * h 
y1 = file_ch01.variables['y_image_bounds'][1] * h 
y2 = file_ch01.variables['y_image_bounds'][0] * h 

# Getting the file time and date
add_seconds = int(file_ch01.variables['time_bounds'][0])
date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
date_formated = date.strftime('%Y-%m-%d %H:%M UTC')
date_file = date.strftime('%Y%m%d%H%M')
year = date.strftime('%Y')
month = date.strftime('%m')
day = date.strftime('%d')
hour = date.strftime('%H')
minutes = date.strftime('%M')

extent = [-80.0, -40.0, -20.0, 10.0]

# Variable to remap
variable = "CMI"

# Call the reprojection funcion
grid = remap(path_ch01, variable, extent, resolution, h, a, b, longitude, x1, y1, x2, y2)
     
# Read the data returned by the function 
data_ch01 = grid.ReadAsArray()
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Read the image
file_ch02 = Dataset(path_ch02)

# Variable to remap
variable = "CMI"

# Call the reprojection funcion
grid = remap(path_ch02, variable, extent, resolution, h, a, b, longitude, x1, y1, x2, y2)
     
# Read the data returned by the function 
data_ch02 = grid.ReadAsArray()
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Read the image
file_ch03 = Dataset(path_ch03)

# Variable to remap
variable = "CMI"

# Call the reprojection funcion
grid = remap(path_ch03, variable, extent, resolution, h, a, b, longitude, x1, y1, x2, y2)
     
# Read the data returned by the function 
data_ch03 = grid.ReadAsArray()

# Calculates the solar zenith angle 
from pyorbital import astronomy
from datetime import datetime

print("Calculating the lons and lats...")
# Create the lats and lons based on the extent
lat = np.linspace(extent[3], extent[1], data_ch01.shape[0])
lon = np.linspace(extent[0], extent[2], data_ch01.shape[1])
xx,yy = np.meshgrid(lon,lat)
lons = xx.reshape(data_ch01.shape[0], data_ch01.shape[1])
lats = yy.reshape(data_ch01.shape[0], data_ch01.shape[1])

# Get the year month day hour and minute to apply the zenith correction
utc_time = datetime(int(year), int(month), int(day), int(hour), int(minutes))
sun_zenith = np.zeros((data_ch01.shape[0], data_ch01.shape[1]))
sun_zenith = astronomy.sun_zenith_angle(utc_time, lons, lats)

# Apply the sun zenith correction
data_ch01 = (data_ch01)/(np.cos(np.deg2rad(sun_zenith)))
data_ch02 = (data_ch02)/(np.cos(np.deg2rad(sun_zenith)))
data_ch03 = (data_ch03)/(np.cos(np.deg2rad(sun_zenith)))

print("Applying the Rayleigh correction...")

# Applying the Rayleigh correction
from pyspectral.rayleigh import Rayleigh     # Atmospherioc correction in the visible spectrum 
from pyorbital.astronomy import get_alt_az
from pyorbital.orbital import get_observer_look

# Satellite height
sat_h = file_ch01.variables['goes_imager_projection'].perspective_point_height

sunalt, suna = get_alt_az(utc_time, lons, lats)
suna = np.rad2deg(suna)
#sata, satel = get_observer_look(sat_lon, sat_lat, sat_alt, vis.attrs['start_time'], lons, lats, 0)
sata, satel = get_observer_look(longitude, 0.0, sat_h, utc_time, lons, lats, 0)
satz = 90 - satel

# Reyleigh Correction
atmosphere = 'us-standard'
aerosol_type = 'rayleigh_only'
rayleigh_key = ('GOES-16','abi', atmosphere, aerosol_type)
corrector = Rayleigh('GOES-16', 'abi', atmosphere=atmosphere, aerosol_type=aerosol_type)

sata = sata % 360.
suna = suna % 360.
ssadiff = np.absolute(suna - sata)
ssadiff = np.minimum(ssadiff, 360 - ssadiff)

red = data_ch02 * 100

refl_cor_band_c01 = corrector.get_reflectance(sun_zenith, satz, ssadiff, 'C01', redband=red)
data_ch01 = data_ch01 - (refl_cor_band_c01 / 100)

refl_cor_band_c02 = corrector.get_reflectance(sun_zenith, satz, ssadiff, 'C02', redband=red)
data_ch02 = data_ch02 - (refl_cor_band_c02 / 100)
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# RGB Components
R = data_ch02
G = (data_ch01 + data_ch02) / 2 * 0.93 + 0.07 * data_ch03 
B = data_ch01

# Apply the CIRA Strech
band_data = R 
log_root = np.log10(0.0223)
denom = (1.0 - log_root) * 0.75
band_data *= 0.01
band_data = band_data.clip(np.finfo(float).eps)
band_data = np.log10(band_data)
band_data -= log_root
band_data /= denom
R  = 1 + band_data
#print (R.shape)

band_data = G
log_root = np.log10(0.0223)
denom = (1.0 - log_root) * 0.75
band_data *= 0.01
band_data = band_data.clip(np.finfo(float).eps)
band_data = np.log10(band_data)
band_data -= log_root
band_data /= denom
G = 1 + band_data
#print (G.shape)

band_data = B
log_root = np.log10(0.0223)
denom = (1.0 - log_root) * 0.75
band_data *= 0.01
band_data = band_data.clip(np.finfo(float).eps)
band_data = np.log10(band_data)
band_data -= log_root
band_data /= denom
B = 1 + band_data
#print (B.shape)

# Create the RGB
RGB = np.stack([R, G, B], axis=2)		
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Product Name
product = "true-color" 
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------    
# Plot configuration
plot_config = {
"resolution": band_resolution_km, 
"dpi": 150, 
"states_color": 'white', "states_width": data_ch01.shape[0] * 0.00018, 
"countries_color": 'white', "countries_width": data_ch01.shape[0] * 0.00100,
"continents_color": 'white', "continents_width": data_ch01.shape[0] * 0.00025,
"grid_color": 'white', "grid_width": data_ch01.shape[0] * 0.00025, "grid_interval": 10.0,
"vmin": 0, "vmax": 1, "cmap": 'jet',
"title_text": "GOES-" + satellite[1:3] + " True color RGB ", "title_size": int(data_ch01.shape[1] * 0.005), "title_x_offset": int(data_ch01.shape[1] * 0.01), "title_y_offset": data_ch01.shape[0] - int(data_ch01.shape[0] * 0.016), 
"thick_interval": 0, "cbar_labelsize": int(data_ch01.shape[0] * 0.005), "cbar_labelpad": -int(data_ch01.shape[0] * 0.0),
"file_name_id_1": satellite,  "file_name_id_2": product
}
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
fig = plt.figure(figsize=(data_ch01.shape[1]/float(plot_config["dpi"]), data_ch01.shape[0]/float(plot_config["dpi"])), dpi=plot_config["dpi"])
  
# Define the projection
proj = ccrs.PlateCarree()

# Use the PlateCarree projection in cartopy
ax = plt.axes([0, 0, 1, 1], projection=proj)
ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Define the image extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]

# Plot the image
img = ax.imshow(RGB, origin='upper', extent=img_extent, zorder=1)

# Add countries
shapefile = list(shpreader.Reader(f'{dir_shapefiles}ne_50m_admin_0_countries.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=plot_config["countries_color"], facecolor='none', linewidth=plot_config["countries_width"], zorder=3)
# Add continents
shapefile = list(shpreader.Reader(f'{dir_shapefiles}ne_10m_coastline.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=plot_config["continents_color"], facecolor='none', linewidth=plot_config["continents_width"], zorder=4)

# Add coastlines, borders and gridlines
ax.gridlines(color=plot_config["grid_color"], alpha=0.5, linestyle='--', linewidth=plot_config["grid_width"], xlocs=np.arange(-180, 180, plot_config["grid_interval"]), ylocs=np.arange(-180, 180, plot_config["grid_interval"]), draw_labels=False, zorder=5)
  
# Add a title
plt.annotate(plot_config["title_text"] + " " + date_formated , xy=(plot_config["title_x_offset"], plot_config["title_y_offset"]), xycoords='figure pixels', fontsize=plot_config["title_size"], fontweight='bold', color='white', bbox=dict(boxstyle="round",fc=(0.0, 0.0, 0.0), ec=(1., 1., 1.)), zorder=6)

# Save the image
plt.savefig(dir_out + plot_config["file_name_id_1"] + "_" + plot_config["file_name_id_2"] + "_" + date_file + '.png', facecolor='black')

# Listar todos os arquivos no diretório
arquivos = os.listdir(dir_in)

# Filtrar os arquivos com extensão .nc.aux.xml
arquivos_aux = [arquivo for arquivo in arquivos if arquivo.endswith('.nc.aux.xml')]

# Excluir os arquivos auxiliares
for arquivo in arquivos_aux:
    caminho_completo = os.path.join(dir_in, arquivo)
    os.remove(caminho_completo)

print('Total processing time:', round((t.time() - start),2), 'seconds.') 
