
#import matplotlib
#matplotlib.use('Agg')
from netCDF4 import Dataset                                  
from mpl_toolkits.axes_grid1.inset_locator import inset_axes                      
import matplotlib.pyplot as plt                              
import numpy as np                                           
import cartopy, cartopy.crs as ccrs                          
from remap import remap                                      
import warnings
warnings.filterwarnings("ignore")
import datetime                                              
from datetime import timedelta  
import time as t                                             
from remap import loadCPT
from utilities import area_para_recorte
from utilities import adicionando_shapefile
from utilities import adicionando_linhas
from utilities import adicionando_descricao_imagem
from utilities import adicionando_logos
from utilities import apply_cira_stretch
from utilities import applying_rayleigh_correction
from utilities import calculating_lons_lats
from dirs import get_dirs
import logging
from osgeo import gdal, osr, ogr 
import os
import warnings
warnings.filterwarnings("ignore")

def process_truecolor(rgb_type, v_extent, ch01, ch02, ch03, ch13):
    global dir_maps
    start = t.time()  
    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)
    # Variable to remap
    variable = "CMI"
    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------

    # Lê a imagem da banda 01
    file_ch01 = Dataset(ch01)
    # Lê a longitude central
    longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin
    # Lê a data do arquivo
    add_seconds = int(file_ch01.variables['time_bounds'][0])
    date = datetime.datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
    date_file = date.strftime('%Y%m%d_%H%M%S')
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')

    #------------------------------------------------------------------------------------------------------#
    #-------------------------------------------Reprojetando----------------------------------------------#
    #------------------------------------------------------------------------------------------------------#
    # reprojetando band 01
    grid = remap(ch01, variable, extent, resolution)
    # Lê o retorno da função
    data_ch01 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # reprojetando band 02
    grid = remap(ch02, variable, extent, resolution)
    # Lê o retorno da função
    data_ch02 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # reprojetando band 03
    grid = remap(ch03, variable, extent, resolution)
    # Lê o retorno da função 
    data_ch03 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
     # reprojetando band 13
    grid = remap(ch13, variable, extent, resolution)
    # Lê o retorno da função
    data_ch13 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    
    # Calculando correção zenith
    utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03 = calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03)

    # Aplicando a correção de Rayleigh
    data_ch01, data_ch02 = applying_rayleigh_correction(file_ch01, utc_time, lons, lats, sun_zenith, data_ch01, data_ch02, longitude)

    # Calculando as cores verdadeiras (True color)
    R = data_ch02
    G = (data_ch01 + data_ch02) / 2 * 0.93 + 0.07 * data_ch03 
    B = data_ch01
    # Aplicando o estiramento CIRA
    R = apply_cira_stretch(R)
    G = apply_cira_stretch(G)
    B = apply_cira_stretch(B)

    # Create the RGB
    RGB = np.stack([R, G, B], axis=2)		
    
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # If zenith angle is greater than 85°, the composite pixel is zero
    RGB[sun_zenith > 85] = 0
    # Create the mask for the regions with zero
    mask = (RGB == [0.0,0.0,0.0]).all(axis=2)
    # Apply the mask to overwrite the pixels
    RGB[mask] = [0,0,0]
    # Create the fading transparency between the regions with the sun zenith angle of 75° and 85°
    alphas = sun_zenith / 100
    min_sun_angle = 0.75
    max_sun_angle = 0.85
    # Normalize the transparency mask
    alphas = ((alphas - max_sun_angle) / (min_sun_angle - max_sun_angle))
    RGB = np.dstack((RGB, alphas))
    
    # if v_extent == 'sp':   #Usar este caso preciso mudar a área
    #     raster = gdal.Open(f'{dir_maps}BlackMarble_2016_B2_geo.tif')
    # else:
    #     raster = gdal.Open(f'{dir_maps}BlackMarble_2016_3km_geo.tif')
    # min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]
    # raster = gdal.Translate(f'{dir_maps}brasil.tif', raster, projWin = [min_lon, max_lat, max_lon, min_lat])

    if v_extent == 'sp':
        raster = gdal.Open(f'{dir_maps}sp_night.tif')
    else:
        raster = gdal.Open(f'{dir_maps}brasil_night.tif')
   
    #lendo o RGB 
    array = raster.ReadAsArray()
    R_night = array[0,:,:].astype(float) / 255
    G_night = array[1,:,:].astype(float) / 255
    B_night = array[2,:,:].astype(float) / 255
    
    R_night[R_night==4] = 0
    G_night[G_night==5] = 0
    B_night[B_night==15] = 0
    
   
    rgb_night = np.stack([R_night, G_night, B_night], axis=2)
    
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  
    
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    
    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-16 Natural True Color + 10.3 µm {date_img}'
    institution = "CEPAGRI - UNICAMP"

    d_p_i = 150
    # Formatando a descricao a ser plotada na imagem
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')
    
    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    ax.imshow(rgb_night, extent=img_extent)
    
    #band13
    data1 = data_ch13
    data1 = np.maximum(data1, 90)
    data1 = np.minimum(data1, 313)
    data1 = (data1-90)/(313-90)
    data1 = 1 - data1

    # Plotando 
    ax.imshow(data1, cmap='gray', vmin=0.1, vmax=0.25, alpha = 0.1, origin='upper', extent=img_extent)

    # Plotando a imagem  # TrueColor
    ax.imshow(RGB, origin='upper', extent=img_extent)
    
    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)

    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)

    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)

    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}{rgb_type}/{rgb_type}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

    print(f'Total processing time: {round((t.time() - start),2)} seconds.')

def iniciar_processo_truecolor(p_br, p_sp, bands, new_bands):
    # Checagem se e possivel gerar imagem TrueColor
    if bands['17']:        
        # Pegando nome dos produtos 01, 02, 03, 13
        ch01 = new_bands['01']
        ch02 = new_bands['02']
        ch03 = new_bands['03']
        ch13 = new_bands['13']
        
        # Pegando nome e local e do produto
        ch01 = f'{dir_in}band01/{ch01}'
        ch02 = f'{dir_in}band02/{ch02}'
        ch03 = f'{dir_in}band03/{ch03}'
        ch13 = f'{dir_in}band13/{ch13}'
        
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS TRUECOLOR WITH NIGHT "BR"...')
            # Tenta realizar o processamento da imagem
            try:
                # Inicia a funcao de processamento
                process_truecolor("truecolor", "br", ch01, ch02, ch03, ch13)
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")

        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_sp:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS TRUECOLOR WITH NIGHT "BR"...')
            # Tenta realizar o processamento da imagem
            try:
                # Inicia a funcao de processamento
                process_truecolor("truecolor", "sp", ch01, ch02, ch03, ch13)
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")


dirs = get_dirs()
# Importando dirs do modulo dirs.py
dir_maps = dirs['dir_maps']
dir_in = dirs['dir_in']
dir_out = dirs['dir_out']
dir_shapefiles = dirs['dir_shapefiles']
dir_colortables = dirs['dir_colortables']
dir_logos = dirs['dir_logos']
dir_out = dirs['dir_out']

bands = {}
bands['17'] = True
p_br = False
p_sp = True
new_bands = { '01': f'OR_ABI-L2-CMIPF-M6C01_G16_s20233050910206_e20233050919514_c20233050919576.nc', 
              '02': f'OR_ABI-L2-CMIPF-M6C02_G16_s20233050910206_e20233050919514_c20233050919567.nc',
              '03': f'OR_ABI-L2-CMIPF-M6C03_G16_s20233050910206_e20233050919514_c20233050919566.nc',
              '13': f'OR_ABI-L2-CMIPF-M6C13_G16_s20233050910206_e20233050919527_c20233050919586.nc'
              }

iniciar_processo_truecolor(p_br, p_sp, bands, new_bands)