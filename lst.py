# Required modules
from netCDF4 import Dataset      # Read / Write NetCDF4 files
import matplotlib.pyplot as plt  # Plotting library
from dirs import get_dirs
from datetime import datetime, timedelta
from utilities import area_para_recorte
from utilities import adicionando_shapefile
from utilities import adicionando_linhas
from utilities import adicionando_descricao_imagem
from utilities import adicionando_logos
import time
import cartopy, cartopy.crs as ccrs    
from remap import remap
import matplotlib
matplotlib.use('Agg')# Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt  # Plotagem de resultados, textos, logos, etc.
import cartopy  # Inserir mapas, shapefiles, paralelos, meridianos, latitudes, longitudes, etc.
import cartopy.crs as ccrs  # Utilitario para sistemas de referência e coordenadas
from osgeo import gdal  # Utilitario para a biblioteca GDAL
from osgeo import osr  # Utilitario para a biblioteca GDAL
from netCDF4 import Dataset  # Utilitario para a biblioteca NetCDF4
import numpy as np  # Suporte para arrays e matrizes multidimensionais, com diversas funções matemáticas para trabalhar com estas estruturas
import time  # Utilitario para trabalhar com tempos
from mpl_toolkits.basemap import Basemap
import cartopy.feature as cfeature # features

dirs = get_dirs()
dir_in = dirs['dir_in']
dir_out = dirs['dir_out']
dir_colortables = dirs['dir_colortables']

# Captura a hora para contagem do tempo de processamento da imagem
processing_start_time = time.time()

file = f"{dir_in}lst/OR_ABI-L2-LST2KMF-M6_G16_s20232991000210_e20232991009518_c20232991010367.nc"

v_extent = 'br'

# Area de interesse para recorte
extent, resolution = area_para_recorte(v_extent)

satellite = '16'
variable = 'LST'

grid = remap(file, variable, extent, resolution)

# Lê o retorno da função
data_lst = grid.ReadAsArray()- 273.15

# Abrindo imagem com a biblioteca GDAL
raw = gdal.Open(f'NETCDF:{file}:' + 'LST', gdal.GA_ReadOnly)
metadata = raw.GetMetadata()
dtime = metadata.get('NC_GLOBAL#time_coverage_start')
date = (datetime.strptime(dtime, '%Y-%m-%dT%H:%M:%S.%fZ'))
date_img = date.strftime('%d-%b-%Y %H:%M UTC')
date_file = date.strftime('%Y%m%d_%H%M%S')

# Define a temperatura minima
min_temp = -26
# Mask values less than -20 degrees and values close to NaN
data_lst = np.ma.masked_where((data_lst < min_temp), data_lst)

# Formatando a descricao a ser plotada na imagem
description = f' GOES-{satellite} Land Surface Temperature {date_img}'
institution = "CEPAGRI - UNICAMP"

# Definindo tamanho da imagem de saida
d_p_i = 150
fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

# Utilizando projecao geoestacionaria no cartopy
ax = plt.axes(projection=ccrs.PlateCarree())

# Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat

# facecolor='DarkSeaGreen'
ax.add_feature(cfeature.LAND, facecolor='SlateGrey', zorder=1)

# facecolor= cornflowerblue deepskyblue dodgerblue
ax.add_feature(cfeature.OCEAN, facecolor='SteelBlue', zorder=2)

# Plotando a imagem
img = ax.imshow(data_lst, origin='upper',vmin=-25, vmax=60,extent=img_extent,  zorder=2, cmap='Spectral_r')

# Adicionando barra da paleta de cores de acordo com o canal
cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[-20, 10, 30, 50])
cb.ax.set_xticklabels(['-20', '10', '30','50'])
cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de coreses

# Adicionando  linhas dos litorais
adicionando_linhas(ax)

# Adicionando descricao da imagem
adicionando_descricao_imagem(description, institution, ax, fig)

# Adicionando os logos
adicionando_logos(fig)

lst = 'lst'
# Salvando a imagem de saida
plt.savefig(f'{dir_out}{lst}/{lst}__{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

# Fecha a janela para limpar a memoria
plt.close()

 