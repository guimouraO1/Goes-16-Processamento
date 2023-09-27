# Required modules
#--------------------------------
#to run in a pure text terminal:
#import matplotlib
#matplotlib.use('Agg')
#--------------------------------
from netCDF4 import Dataset                                  # Read / Write NetCDF4 files
from mpl_toolkits.axes_grid1.inset_locator import inset_axes # Add a child inset axes to this existing axes.
from datetime import datetime, timedelta                     # Library to convert julian day to dd-mm-yyyy                            # Import the CPT convert function
from matplotlib.colors import LinearSegmentedColormap        # Linear interpolation for color maps
import matplotlib.pyplot as plt                              # Plotting library
import numpy as np                                           # Scientific computing with Python
import cartopy, cartopy.crs as ccrs                          # Plot maps
import cartopy.io.shapereader as shpreader                   # Import shapefiles
import time as t                                             # Time access and conversion
import math                                                  # Import math                                                # Unix style pathname pattern expansion
import sys                                                   # Import the "system specific parameters and functions" module
import os 												     # Miscellaneous operating system interfaces
from os.path import dirname, abspath                         # Return a normalized absolutized version of the pathname path                                            # To check which OS is being used
from pyorbital import astronomy                      # Import GDAL                            # Update the HTML animation 
from remap import remap                                      # Import the Remap function 
import warnings
warnings.filterwarnings("ignore")


import colorsys							 # Conversions between color systems

def loadCPT(path):

    try:
        f = open(path)
    except:
        print ("File ", path, "not found")
        return None

    lines = f.readlines()

    f.close()

    x = np.array([])
    r = np.array([])
    g = np.array([])
    b = np.array([])

    colorModel = 'RGB'

    for l in lines:
        ls = l.split()
        if l[0] == '#':
            if ls[-1] == 'HSV':
                colorModel = 'HSV'
                continue
            else:
                continue
        if ls[0] == 'B' or ls[0] == 'F' or ls[0] == 'N':
            pass
        else:
            x=np.append(x,float(ls[0]))
            r=np.append(r,float(ls[1]))
            g=np.append(g,float(ls[2]))
            b=np.append(b,float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

        x=np.append(x,xtemp)
        r=np.append(r,rtemp)
        g=np.append(g,gtemp)
        b=np.append(b,btemp)

    if colorModel == 'HSV':
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
        r[i] = rr ; g[i] = gg ; b[i] = bb

    if colorModel == 'RGB':
        r = r/255.0
        g = g/255.0
        b = b/255.0

    xNorm = (x - x[0])/(x[-1] - x[0])

    red   = []
    blue  = []
    green = []

    for i in range(len(x)):
        red.append([xNorm[i],r[i],r[i]])
        green.append([xNorm[i],g[i],g[i]])
        blue.append([xNorm[i],b[i],b[i]])

    colorDict = {'red': red, 'green': green, 'blue': blue}
    #print(red)

    return colorDict

# Start the time counter
print('Script started.')
start = t.time()  

# Band 02 path
path_ch02 = '/home/guimoura/Documentos/Goes-16-Processamento/goes/band02/OR_ABI-L2-CMIPF-M6C02_G16_s20232701130209_e20232701139517_c20232701139567.nc'
path_ch13 = '/home/guimoura/Documentos/Goes-16-Processamento/goes/band13/OR_ABI-L2-CMIPF-M6C13_G16_s20232701130209_e20232701139529_c20232701139586.nc'

# Read the image
file_ch02 = Dataset(path_ch02)

# Read the satellite 
satellite = getattr(file_ch02, 'platform_ID')

# Desired resolution
resolution = 4

# Read the central longitude
longitude = file_ch02.variables['goes_imager_projection'].longitude_of_projection_origin


# Getting the file time and date
add_seconds = int(file_ch02.variables['time_bounds'][0])
date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
date_formated = date.strftime('%Y-%m-%d %H:%M UTC')
date_file = date.strftime('%Y%m%d%H%M')
year = date.strftime('%Y')
month = date.strftime('%m')
day = date.strftime('%d')
hour = date.strftime('%H')
minutes = date.strftime('%M')

# Choose the visualization extent (min lon, min lat, max lon, max lat)
extent = [-90.0, -40.0, -20.0, 10.0]

# Variable to remap
variable = "CMI"

# Call the reprojection funcion
grid = remap(path_ch02, variable, extent, resolution)
     
# Read the data returned by the function 
data_ch02 = grid.ReadAsArray()

# Load the Channel 02 contrast curve
contrast_curve = [0.00000, 0.02576, 0.05148, 0.07712, 0.10264, 0.12799, 0.15313, 0.17803, 0.20264, 0.22692, 0.25083, 0.27432, 0.29737, 0.31991, 0.34193, 0.36336,
0.38418, 0.40433, 0.42379, 0.44250, 0.46043, 0.47754, 0.49378, 0.50911, 0.52350, 0.53690, 0.54926, 0.56055, 0.57073, 0.57976, 0.58984, 0.59659,
0.60321, 0.60969, 0.61604, 0.62226, 0.62835, 0.63432, 0.64016, 0.64588, 0.65147, 0.65694, 0.66230, 0.66754, 0.67267, 0.67768, 0.68258, 0.68738,
0.69206, 0.69664, 0.70112, 0.70549, 0.70976, 0.71394, 0.71802, 0.72200, 0.72589, 0.72968, 0.73339, 0.73701, 0.74055, 0.74399, 0.74736, 0.75065,
0.75385, 0.75698, 0.76003, 0.76301, 0.76592, 0.76875, 0.77152, 0.77422, 0.77686, 0.77943, 0.78194, 0.78439, 0.78679, 0.78912, 0.79140, 0.79363,
0.79581, 0.79794, 0.80002, 0.80206, 0.80405, 0.80600, 0.80791, 0.80978, 0.81162, 0.81342, 0.81518, 0.81692, 0.81862, 0.82030, 0.82195, 0.82358,
0.82518, 0.82676, 0.82833, 0.82987, 0.83140, 0.83292, 0.83442, 0.83592, 0.83740, 0.83888, 0.84036, 0.84183, 0.84329, 0.84476, 0.84623, 0.84771,
0.84919, 0.85068, 0.85217, 0.85368, 0.85520, 0.85674, 0.85829, 0.85986, 0.86145, 0.86306, 0.86469, 0.86635, 0.86803, 0.86974, 0.87149, 0.87326,
0.87500, 0.87681, 0.87861, 0.88038, 0.88214, 0.88388, 0.88560, 0.88730, 0.88898, 0.89064, 0.89228, 0.89391, 0.89552, 0.89711, 0.89868, 0.90023,
0.90177, 0.90329, 0.90479, 0.90627, 0.90774, 0.90919, 0.91063, 0.91205, 0.91345, 0.91483, 0.91620, 0.91756, 0.91890, 0.92022, 0.92153, 0.92282,
0.92410, 0.92536, 0.92661, 0.92784, 0.92906, 0.93027, 0.93146, 0.93263, 0.93380, 0.93495, 0.93608, 0.93720, 0.93831, 0.93941, 0.94050, 0.94157,
0.94263, 0.94367, 0.94471, 0.94573, 0.94674, 0.94774, 0.94872, 0.94970, 0.95066, 0.95162, 0.95256, 0.95349, 0.95441, 0.95532, 0.95622, 0.95711,
0.95799, 0.95886, 0.95973, 0.96058, 0.96142, 0.96225, 0.96307, 0.96389, 0.96469, 0.96549, 0.96628, 0.96706, 0.96783, 0.96860, 0.96936, 0.97010,
0.97085, 0.97158, 0.97231, 0.97303, 0.97374, 0.97445, 0.97515, 0.97584, 0.97653, 0.97721, 0.97789, 0.97856, 0.97922, 0.97988, 0.98053, 0.98118,
0.98182, 0.98246, 0.98309, 0.98372, 0.98435, 0.98497, 0.98559, 0.98620, 0.98681, 0.98741, 0.98802, 0.98862, 0.98921, 0.98980, 0.99039, 0.99098,
0.99157, 0.99215, 0.99273, 0.99331, 0.99389, 0.99446, 0.99503, 0.99561, 0.99618, 0.99675, 0.99732, 0.99788, 0.99845, 0.99902, 0.99959, 1.000000]
 
# Convert the contrast curve to a numpy array
curve = np.array(contrast_curve)

# Minimum and maximum reflectances
VISmin = 0.0
VISmax = 1.0
# Anything that is below the min or greater than max, keep min and max
data_ch02[data_ch02 > VISmax] = VISmax
# Convert to 0 - 255
data_ch02 = ((data_ch02 - VISmin) / (VISmax - VISmin)) * 255
# Convert to int
data_ch02 = data_ch02.astype(int)
# Apply the contrast curve
data_ch02 = curve[data_ch02] * 255
# Convert to int
data_ch02 = data_ch02.astype(int)

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

# Read the image
file_ch13 = Dataset(path_ch13)

# Read the resolution
band_resolution_km = getattr(file_ch13, 'spatial_resolution')
band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

# Division factor to reduce image size
f = math.ceil(float(resolution / band_resolution_km))

# Variable to remap
variable = "CMI"

# Call the reprojection funcion
grid = remap(path_ch13, variable, extent, resolution)
     
# Read the data returned by the function 
data_ch13 = grid.ReadAsArray()
data_ch13_original = data_ch13

# Minimum and maximum reflectances
IRmin = 89.62
IRmax = 341.27
# Anything that is below the min or greater than max, keep min and max
data_ch13[data_ch13 > IRmax] = IRmax
# Convert to 0 - 255
data_ch13 = ((data_ch13 - IRmin) / (IRmax - IRmin)) * 255
# Convert to int
data_ch13 = data_ch13.astype(int)

#print("Reading the False Color Look Up Table...")
 
import matplotlib.image as mpimg
# Open the False Color Look Up Table
img = mpimg.imread('/home/guimoura/Documentos/Goes-16-Processamento/colortables/wx-star.com_GOES-R_ABI_False-Color-LUT_sat20.png')
# Flip the array (for some reason, is necessary to flip the LUT horizontally)
img = np.fliplr(img)
# Apply the look up table based on the Band 02 reflectances and Band 13 Brightness Temperatures
data = img[data_ch02,data_ch13,0:3]

# Eliminate values outside the globe
#mask = (data == [data[0,0]])
#data[mask] = np.nan


#print(data_ch02.shape[0])
lons = np.arange(extent[0], extent[2], ((extent[2] - extent[0]) / data_ch02.shape[1]))
print(lons)
lats = np.arange(extent[1], extent[3], ((extent[3] - extent[1]) / data_ch02.shape[0]))
print(lats)

lons, lats = np.meshgrid(lons, lats)
lats = np.flipud(lats)

#print("Calculating the sun zenith angle...")
utc_time = datetime(int(year), int(month), int(day), int(hour), int(minutes))
sun_zenith = np.zeros((data_ch02.shape[0], data_ch02.shape[1]))
sun_zenith = astronomy.sun_zenith_angle(utc_time, lons[:,:][::1,::1], lats[:,:][::1,::1])


data[sun_zenith > 85] = [0.0,0.0,0.0]
mask = (data == [0.0,0.0,0.0]).all(axis=2)
#apply the mask to overwrite the pixels
data[mask] = [0.0,0.0,0.0]
# Create the fading transparency between the regions with the
# sun zenith angle of 75° and 85°
alphas = sun_zenith / 100
min_sun_angle = 0.75
max_sun_angle = 0.85
alphas = ((alphas - max_sun_angle) / (min_sun_angle - max_sun_angle))
#alpha = ~np.all(data == [0.0,0.0,0.0], axis=2)
data = np.dstack((data, alphas))

# Choose the visualization extent (min lon, min lat, max lon, max lat)
extent = [-90.0, -40.0, -20.0, 10.0]

extent_night = [-156.29, -81.32, 6.29, 81.32]

extent_night[0] = extent[0] - 10
extent_night[1] = extent[1] - 10
extent_night[2] = extent[2] + 10
extent_night[3] = extent[3] + 10

# Plot configuration
plot_config = {
"resolution": band_resolution_km, 
"dpi": 150, 
"states_color": 'white', "states_width": data_ch02.shape[0] * 0.00006, 
"countries_color": 'gold', "countries_width": data_ch02.shape[0] * 0.00012,
"continents_color": 'gold', "continents_width": data_ch02.shape[0] * 0.00025,
"grid_color": 'white', "grid_width": data_ch02.shape[0] * 0.00025, "grid_interval": 10.0,
"title_text": "GOES-" + satellite[1:3] + " FALSE COLOR ", "title_size": int(data_ch02.shape[1] * 0.005), "title_x_offset": int(data_ch02.shape[1] * 0.01), "title_y_offset": data_ch02.shape[0] - int(data_ch02.shape[0] * 0.016), 
"file_name_id_1": satellite,  "file_name_id_2": 'product'
}

# Choose the plot size (width x height, in inches)
fig = plt.figure(figsize=(data_ch02.shape[1]/float(plot_config["dpi"]), data_ch02.shape[0]/float(plot_config["dpi"])), dpi=plot_config["dpi"])
  
# Define the projection
proj = ccrs.PlateCarree()

# Use the PlateCarree projection in cartopy
ax = plt.axes([0, 0, 1, 1], projection=proj)
ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Define the image extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]

#print("First layer...")
# Apply range limits for clean IR channel
data1 = data_ch13_original
data1 = np.maximum(data1, 90)
data1 = np.minimum(data1, 313)
# Normalize the channel between a range
data1 = (data1-90)/(313-90)
# Invert colors
data1 = 1 - data1
img2 = ax.imshow(data1, cmap='gray', vmin=0.1, vmax=0.25, alpha = 0.6, origin='upper', extent=img_extent, zorder=2)

#print("Second layer...")
# SECOND LAYER
data2 = data1
data2[data2 < 0.20] = np.nan
img3 = ax.imshow(data2, cmap='gray', vmin=0.15, vmax=0.30, alpha = 1.0, origin='upper', extent=img_extent, zorder=3)

# Converts a CPT file to be used in Python
cpt = loadCPT('/home/guimoura/Documentos/Goes-16-Processamento/colortables/IR4AVHRR6.cpt')   
cmap = LinearSegmentedColormap('cpt', cpt) 
#print("Third layer...")
# THIRD LAYER
data3 = data_ch13_original
data3 = data_ch13_original - 273.15
data3[np.logical_or(data3 < -80, data3 > -28)] = np.nan
img4 = ax.imshow(data3, cmap=cmap, vmin=-103, vmax=84, alpha=1.0, origin='upper', extent=img_extent, zorder=4)

# Plot the image
ax.imshow(data, origin='upper', extent=img_extent, zorder=5)

# Save the image
plt.savefig('/home/guimoura/Documentos/Goes-16-Processamento/' + plot_config["file_name_id_1"] + "_" + plot_config["file_name_id_2"] + "_" + date_file + '.png', facecolor='black')#, bbox_inches='tight', pad_inches=0, facecolor='black')



print('Total processing time:', round((t.time() - start),2), 'seconds.') 