from netCDF4 import Dataset
import numpy as np
from osgeo import osr
from osgeo import gdal
import time as t
import warnings
warnings.filterwarnings("ignore")
gdal.PushErrorHandler('CPLQuietErrorHandler') 
import colorsys		
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
def remap(path, variable, extent, resolution):
    
    # Default scale    
    scale = 1
    
    # Default offset
    offset = 0
    
    # GOES Extent (satellite projection) [llx, lly, urx, ury]
    GOES_EXTENT = [-5434894.885056, -5434894.885056, 5434894.885056, 5434894.885056]
    
    # Setup NetCDF driver
    gdal.SetConfigOption('GDAL_NETCDF_BOTTOMUP', 'NO')
        
    if not (variable == "DQF"):              
        # Read scale/offset from file
        scale, offset = getScaleOffset(path, variable) 
      
    connectionInfo = f'NETCDF:\"' + path + '\"://' + variable
	
    # Lendo os metadados do cabecalho
    raw = gdal.Open(connectionInfo, gdal.GA_ReadOnly)          


    # # GOES Spatial Reference System
    # sourcePrj = osr.SpatialReference()
    # sourcePrj.ImportFromProj4('+proj=geos +h=' + str(h) + ' ' + '+a=' + str(a) + ' ' + '+b=' + str(b) + ' ' + '+lon_0=' + str(longitude) + ' ' + '+sweep=x')

    # # Lat/lon WSG84 Spatial Reference System
    # targetPrj = osr.SpatialReference()
    # targetPrj.ImportFromProj4('+proj=latlong +datum=WGS84')
    
    # GOES-16 Spatial Reference System
    sourcePrj = osr.SpatialReference()
    sourcePrj.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs')
    # Lat/lon WSG84 Spatial Reference System
    targetPrj = osr.SpatialReference()
    targetPrj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Setup projection and geo-transformation
    raw.SetProjection(sourcePrj.ExportToWkt())
    raw.SetGeoTransform(getGeoT(GOES_EXTENT, raw.RasterYSize, raw.RasterXSize))  

    # Compute grid dimension
    KM_PER_DEGREE = 111.32
    sizex = int(((extent[2] - extent[0]) * KM_PER_DEGREE) / resolution)
    sizey = int(((extent[3] - extent[1]) * KM_PER_DEGREE) / resolution)
    
    # Get memory driver
    memDriver = gdal.GetDriverByName('MEM')
   
    # Create grid
    grid = memDriver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)
        
    # Setup projection and geo-transformation
    grid.SetProjection(targetPrj.ExportToWkt())
    grid.SetGeoTransform(getGeoT(extent, grid.RasterYSize, grid.RasterXSize))
    
    
    gdal.ReprojectImage(raw, grid, sourcePrj.ExportToWkt(), targetPrj.ExportToWkt(), gdal.GRA_NearestNeighbour, options=['NUM_THREADS=ALL_CPUS']) 
    

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
    del raw
	
    return grid


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
