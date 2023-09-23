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
from truecolor import adicionando_shapefile
from truecolor import adicionando_linhas
from truecolor import adicionando_descricao_imagem
from truecolor import adicionando_logos
from remap import remap
import time as time    
import logging

###########################################################################
#              Script de Processamento para True Color Goes-16            #
###########################################################################
#  Metodo: De 10 em 10 minutos                                            #
#  Descricao: Processa imagens netCDF4 para criar imagem True Color       #
#  Autor: Guilherme de Moura Oliveira  <guimoura@unicamp.br>              #
#  Data: 22/09/2023                                                       #
#  Atualizacao: 22/09/2023                                                #
###########################################################################

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


# Lê o identificador do satélite
satellite = getattr(file_ch08, 'platform_ID')

variable = "CMI"
v_extent = 'br'
extent, resolution = area_para_recorte(v_extent)

# Read the resolution
band_resolution_km = getattr(file_ch08, 'spatial_resolution')
band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

# Getting the file time and date
add_seconds = int(file_ch08.variables['time_bounds'][0])
date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
date_formated = date.strftime('%Y-%m-%d %H:%M UTC')
date_file = date.strftime('%Y%m%d_%H%M%S')
date_img = date.strftime('%d-%b-%Y %H:%M UTC')
year = date.strftime('%Y')
month = date.strftime('%m')
day = date.strftime('%d')
hour = date.strftime('%H')
minutes = date.strftime('%M')

#------------------------------------------------------------------------------------------------------#
#-------------------------------------------Reprojetando----------------------------------------------#
#------------------------------------------------------------------------------------------------------#
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

# Formatando a descricao a ser plotada na imagem
description = f' GOES-{satellite} Air Mass {date_img}'
institution = "CEPAGRI - UNICAMP"

d_p_i = 150
fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

# Utilizando projecao geoestacionaria no cartopy
ax = plt.axes(projection=ccrs.PlateCarree())

# Adicionando o shapefile dos estados brasileiros
adicionando_shapefile(v_extent, ax)

# Adicionando  linhas dos litorais
adicionando_linhas(ax)

# Define the image extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]

# Plot the image
img = ax.imshow(RGB, origin='upper', extent=img_extent)

# Adicionando descricao da imagem
adicionando_descricao_imagem(description, institution, ax, fig)

# Adicionando os logos
adicionando_logos(fig)

dir_main =  f'/mnt/e/truecolor/' 
dir_out = f'{dir_main}output/'
rgb_type = 'airmass'
# Salvando a imagem de saida
plt.savefig(f'{dir_out}{rgb_type}/{rgb_type}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

# Fecha a janela para limpar a memoria
plt.close()

# Tempo de processamento True color
logging.info('Total processing time Airmass:', round((time.time() - start),2), 'seconds.') 

