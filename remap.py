from netCDF4 import Dataset                  # Read / Write NetCDF4 files
import matplotlib.pyplot as plt              # Plotting library
import cartopy, cartopy.crs as ccrs          # Plot maps
import numpy as np                           # Numpy library
from pyorbital import astronomy              # Compute orbital parameters for satellites from TLE files
from datetime import datetime, timedelta     # Library to convert julian day to dd-mm-yyyy
import math                                  # Import math
from pyproj import Proj                      # Python interface to PROJ (cartographic projections and coordinate transformations library)
from pyspectral.rayleigh import Rayleigh     # Atmospherioc correction in the visible spectrum 
import cartopy.io.shapereader as shpreader   # Import shapefiles
import numpy as np
import time  # Utilitario para trabalhar com tempos
# Calculates the solar zenith angle 
from pyorbital import astronomy
from datetime import datetime

# Suprimir todos os avisos do NumPy
np.seterr(all='ignore')

dir_in = f'/mnt/e/truecolor/'
dir_shapefiles = f'{dir_in}shapefiles/'
dir_colortables = f'{dir_in}colortables/'
dir_logos = f'{dir_in}logos/'
dir_out = f'{dir_in}output/'
# File to read
image1 = f'{dir_in}reprojetadas/CG_ABI-L2-CMIPF-M6C01_G16_s20232631130206_e20232631139514_c20232631139581_reproj_br.nc'
image2 = f'{dir_in}reprojetadas/CG_ABI-L2-CMIPF-M6C02_G16_s20232631130206_e20232631139514_c20232631139575_reproj_br.nc'
image3 = f'{dir_in}reprojetadas/CG_ABI-L2-CMIPF-M6C03_G16_s20232631130206_e20232631139514_c20232631139581_reproj_br.nc'

# Read the file using the NetCDF library
file1 = Dataset(image1)
file2 = Dataset(image2)
file3 = Dataset(image3)

#######################################################################################################
# LICENSE
# Copyright (C) 2021 - INPE - NATIONAL INSTITUTE FOR SPACE RESEARCH - BRAZIL
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU 
# General Public License as published by the Free Software Foundation, either version 3 of the License, 
# or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
# Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. 
# If not, see http://www.gnu.org/licenses/.
#######################################################################################################
__author__ = "Diego Souza"
__copyright__ = "Copyright (C) 2021 - INPE - NATIONAL INSTITUTE FOR SPACE RESEARCH - BRAZIL"
__credits__ = ["Diego Souza", "Douglas Uba"]
__license__ = "GPL"
__version__ = "2.5.0"
__maintainer__ = "Diego Souza"
__email__ = "diego.souza@inpe.br"
__status__ = "Production"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Required modules
from netCDF4 import Dataset
import numpy as np
from osgeo import osr
from osgeo import gdal
import time as t
import warnings
warnings.filterwarnings("ignore")
gdal.PushErrorHandler('CPLQuietErrorHandler') 
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
def exportImage(image,path):
    driver = gdal.GetDriverByName('netCDF')
    return driver.CreateCopy(path,image,0)
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
def getGeoT(extent, nlines, ncols):
    # Compute resolution based on data dimension
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3] , 0, -resy]
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
def getScaleOffset(path, variable):
    nc = Dataset(path, mode='r')
    
    if (variable == "BCM") or (variable == "Phase") or (variable == "Smoke") or (variable == "Dust") or (variable == "Mask") or (variable == "Power"): 
        scale  = 1
        offset = 0     
    else:
        scale = nc.variables[variable].scale_factor
        offset = nc.variables[variable].add_offset
    nc.close()
        
    return scale, offset
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------    
def remap(path, variable, extent, resolution, h, a, b, longitude, x1, y1, x2, y2):
    
    # Default scale    
    scale = 1
    
    # Default offset
    offset = 0
    
    # GOES Extent (satellite projection) [llx, lly, urx, ury]
    GOES_EXTENT = [x1, y1, x2, y2]
    
    # Setup NetCDF driver
    gdal.SetConfigOption('GDAL_NETCDF_BOTTOMUP', 'NO')
        
    if not (variable == "DQF"):              
        # Read scale/offset from file
        scale, offset = getScaleOffset(path, variable) 
      
    connectionInfo = 'HDF5:\"' + path + '\"://' + variable
    
    print(connectionInfo)
	
    # Read the datasat
    raw = gdal.Open(connectionInfo)          
    
    # Define KM_PER_DEGREE
    KM_PER_DEGREE = 111.32

    # GOES Spatial Reference System
    sourcePrj = osr.SpatialReference()
    sourcePrj.ImportFromProj4('+proj=geos +h=' + str(h) + ' ' + '+a=' + str(a) + ' ' + '+b=' + str(b) + ' ' + '+lon_0=' + str(longitude) + ' ' + '+sweep=x')

    # Lat/lon WSG84 Spatial Reference System
    targetPrj = osr.SpatialReference()
    targetPrj.ImportFromProj4('+proj=latlong +datum=WGS84')

    # Setup projection and geo-transformation
    raw.SetProjection(sourcePrj.ExportToWkt())
    raw.SetGeoTransform(getGeoT(GOES_EXTENT, raw.RasterYSize, raw.RasterXSize))  

    # Compute grid dimension
    sizex = int(((extent[2] - extent[0]) * KM_PER_DEGREE) / resolution)
    sizey = int(((extent[3] - extent[1]) * KM_PER_DEGREE) / resolution)
    
    # Get memory driver
    memDriver = gdal.GetDriverByName('MEM')
   
    # Create grid
    grid = memDriver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)
        
    # Setup projection and geo-transformation
    grid.SetProjection(targetPrj.ExportToWkt())
    grid.SetGeoTransform(getGeoT(extent, grid.RasterYSize, grid.RasterXSize))

    # Perform the projection/resampling 
    print ('Remapping...')#, path)
        
    start = t.time()
    
    gdal.ReprojectImage(raw, grid, sourcePrj.ExportToWkt(), targetPrj.ExportToWkt(), gdal.GRA_NearestNeighbour, options=['NUM_THREADS=ALL_CPUS']) 
    
    print ('Remap finished! Time:', round((t.time() - start),2), 'seconds')
               
    # Read grid data
    array = grid.ReadAsArray()
    
    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)
    
    # Read as uint16
    array = array.astype(np.uint16)  
       
    # Apply scale and offset
    array = array * scale + offset

    # Get the raster 
    grid.GetRasterBand(1).WriteArray(array)

	# Close file
    raw = None
	
    return grid
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------