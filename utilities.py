import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset                                  # Lê / Escrever arquivos NetCDF4
import datetime                                              # Biblioteca para trabalhar com datas
from datetime import timedelta, datetime                     # Biblioteca para converter dia juliano para dd-mm-yyyy
import matplotlib.pyplot as plt                              # Biblioteca para plotar gráficos
import matplotlib.colors                                     # Cores do Matplotlib
import numpy as np                                           # Computação científica com Python
import cartopy, cartopy.crs as ccrs                          # Plotar mapas
import cartopy.io.shapereader as shpreader                   # Importar shapefiles     
from pyorbital import astronomy
from pyspectral.rayleigh import Rayleigh                     # Correção atmosférica no espectro visível 
from pyorbital.astronomy import get_alt_az
from pyorbital.orbital import get_observer_look
from dirs import get_dirs
import os
import boto3  # Amazon Web Services (AWS) SDK for Python
from botocore import UNSIGNED  # boto3 config
from botocore.config import Config  # boto3 config
from osgeo import osr
from osgeo import gdal
import warnings
warnings.filterwarnings("ignore")
gdal.PushErrorHandler('CPLQuietErrorHandler') 
import colorsys

dirs = get_dirs()
# Importando dirs do modulo dirs.py
dir_main = dirs['dir_main']
dir_maps = dirs['dir_maps']
dir_in = dirs['dir_in']
dir_out = dirs['dir_out']
dir_shapefiles = dirs['dir_shapefiles']
dir_colortables = dirs['dir_colortables']
dir_logos = dirs['dir_logos']
dir_out = dirs['dir_out']


def adicionando_linhas(ax):
    # Adicionando linhas da costa
    ax.coastlines(resolution='10m', color='cyan', linewidth=0.5, zorder=5)
    # Adicionando linhas das fronteiras
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='cyan', linewidth=0.5, zorder=5)
    # Adicionando paralelos e meridianos
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=0.7, linestyle='--', linewidth=0.2, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5))
    gl.top_labels = False
    gl.right_labels = False


def area_para_recorte(v_extent):
    
    # Area de interesse para recorte
    if v_extent == 'br':
        # Brasil
        extent = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat
    # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 4.0
    elif v_extent == 'sp':
        # São Paulo
        extent = [-53.25, -26.0, -44.0, -19.5]  # Min lon, Min lat, Max lon, Max lat
    # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 1.0
    else:
        extent = [-115.98, -55.98, -25.01, 34.98]  # Min lon, Min lat, Max lon, Max lat
        resolution = 2.0
    return extent, resolution


def adicionando_shapefile(v_extent, ax):
    if v_extent == 'br':
        # Adicionando o shapefile dos estados brasileiros
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='cyan', facecolor='none', linewidth=0.7, zorder=5)
    elif v_extent == 'sp':
        # Adicionando o shapefile dos estados brasileiros e da cidade de Campinas
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='cyan', facecolor='none', linewidth=0.7, zorder=5)
        shapefile = list(shpreader.Reader(dir_shapefiles + 'campinas/campinas.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='cyan', facecolor='none', linewidth=1, zorder=5)


def adicionando_descricao_imagem(description, institution, ax, fig, cruz=False):
    cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
    cax1.patch.set_color('black')  # Alterando a cor do novo eixo
    cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
    if cruz:
        cruzStr = '+'
        cax1.text(0.190, 0.13, cruzStr, color='red', size=12)  # Adicionando símbolo "+"
    cax1.text(0.85, 0.13, institution, color='yellow', size=10)  # Adicionando texto
    cax1.xaxis.set_visible(False)  # Removendo rótulos do eixo X
    cax1.yaxis.set_visible(False)  # Removendo rótulos do eixo Y


def calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03):
        
    year = date.strftime('%Y')
    month = date.strftime('%m')
    day = date.strftime('%d')
    hour = date.strftime('%H')
    minutes = date.strftime('%M')
    
    # Criar as latitudes e longitudes com base na extensão
    lat = np.linspace(extent[3], extent[1], data_ch01.shape[0])
    lon = np.linspace(extent[0], extent[2], data_ch01.shape[1])
    xx,yy = np.meshgrid(lon,lat)
    lons = xx.reshape(data_ch01.shape[0], data_ch01.shape[1])
    lats = yy.reshape(data_ch01.shape[0], data_ch01.shape[1])

    # Obter o ano, mês, dia, hora e minuto para aplicar a correção zenital
    utc_time = datetime(int(year), int(month), int(day), int(hour), int(minutes))
    sun_zenith = np.zeros((data_ch01.shape[0], data_ch01.shape[1]))
    sun_zenith = astronomy.sun_zenith_angle(utc_time, lons, lats)

    # Aplicar a correção zenital do sol
    data_ch01 = (data_ch01)/(np.cos(np.deg2rad(sun_zenith)))
    data_ch02 = (data_ch02)/(np.cos(np.deg2rad(sun_zenith)))
    data_ch03 = (data_ch03)/(np.cos(np.deg2rad(sun_zenith)))

    return utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03


def applying_rayleigh_correction(file_ch01, utc_time, lons, lats, sun_zenith, data_ch01, data_ch02, longitude):
    # Altitude do satélite
    sat_h = file_ch01.variables['goes_imager_projection'].perspective_point_height

    sunalt, suna = get_alt_az(utc_time, lons, lats)
    suna = np.rad2deg(suna)
    #sata, satel = get_observer_look(sat_lon, sat_lat, sat_alt, vis.attrs['start_time'], lons, lats, 0)
    sata, satel = get_observer_look(longitude, 0.0, sat_h, utc_time, lons, lats, 0)
    satz = 90 - satel

    # Correção de Rayleigh
    atmosphere = 'us-standard'
    aerosol_type = 'rayleigh_only'
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
    
    return data_ch01, data_ch02 


def adicionando_logos(fig):
    # Adicionando os logotipos
    logo_noaa = plt.imread(dir_logos + 'NOAA_Logo.png')  # Lendo o arquivo do logotipo
    logo_goes = plt.imread(dir_logos + 'GOES_Logo.png')  # Lendo o arquivo do logotipo
    logo_cepagri = plt.imread(dir_logos + 'CEPAGRI-Logo.png')  # Lendo o arquivo do logotipo
    fig.figimage(logo_noaa, 32, 240, zorder=3, alpha=0.6, origin='upper')  # Plotando logotipo
    fig.figimage(logo_goes, 10, 160, zorder=3, alpha=0.6, origin='upper')  # Plotando logotipo
    fig.figimage(logo_cepagri, 10, 80, zorder=3, alpha=0.8, origin='upper')  # Plotando logotipo


def apply_cira_stretch(band_data):
    
    log_root = np.log10(0.0223)
    denom = (1.0 - log_root) * 0.75
    band_data *= 0.01
    band_data = band_data.clip(np.finfo(float).eps)
    band_data = np.log10(band_data)
    band_data -= log_root
    band_data /= denom
    return 1 + band_data


def download_prod(yyyymmddhhmn, product_name, path_dest):
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6_G16_s{year}{day_of_year}{hour}{min}'
    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmn}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}{file_name}.nc exists')
            else:
                # print(f'Downloading file {path_dest}{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}'


def getGeoT(extent, nlines, ncols):
    # Compute resolution based on data dimension
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3] , 0, -resy]


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
