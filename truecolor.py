import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset                                  # Read / Write NetCDF4 files
from datetime import datetime, timedelta                     # Library to convert julian day to dd-mm-yyyy
import matplotlib.pyplot as plt                              # Plotting library
import matplotlib.colors                                     # Matplotlib colors
import numpy as np                                           # Scientific computing with Python
import cartopy, cartopy.crs as ccrs                          # Plot maps
import cartopy.io.shapereader as shpreader                   # Import shapefiles
import time as t                                             # Time access and conversion
import math                                                  # Import math
import os 												     # Miscellaneous operating system interfaces
from os.path import dirname, abspath                         # Return a normalized absolutized version of the pathname path 
from remap import remap                                      # Import the Remap function
import warnings
from pyorbital import astronomy
from datetime import datetime
from pyspectral.rayleigh import Rayleigh     # Atmospherioc correction in the visible spectrum 
from pyorbital.astronomy import get_alt_az
from pyorbital.orbital import get_observer_look

warnings.filterwarnings("ignore")


def adicionando_linhas(ax):
    # Adicionando  linhas dos litorais
    ax.coastlines(resolution='10m', color='cyan', linewidth=0.5)
    # Adicionando  linhas das fronteiras
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='cyan', linewidth=0.5)
    # Adicionando  paralelos e meridianos
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=0.7, linestyle='--', linewidth=0.2, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5))
    gl.top_labels = False
    gl.right_labels = False


def adicionando_shapefile(v_extent, ax):
    if v_extent == 'br':
        # Adicionando o shapefile dos estados brasileiros
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='cyan', facecolor='none', linewidth=0.7)
    elif v_extent == 'sp':
        # Adicionando o shapefile dos estados brasileiros e cidade de campinas
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='cyan', facecolor='none', linewidth=0.7)
        shapefile = list(shpreader.Reader(dir_shapefiles + 'campinas/campinas').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='yellow', facecolor='none', linewidth=1)


def adicionando_descricao_imagem(description, institution, ax, fig, cruz=False):
    # Criando novos eixos de acordo com a posicao da imagem
    cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
    cax1.patch.set_color('black')  # Alterando a cor do novo eixo
    cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
    if cruz:
        cruzStr = '+'
        cax1.text(0.190, 0.13, cruzStr, color='red', size=12)  # Adicionando simbolo "+"
    cax1.text(0.85, 0.13, institution, color='yellow', size=10)  # Adicionando texto
    cax1.xaxis.set_visible(False)  # Removendo rotulos do eixo X
    cax1.yaxis.set_visible(False)  # Removendo rotulos do eixo Y


def calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03):
        
    year = date.strftime('%Y')
    month = date.strftime('%m')
    day = date.strftime('%d')
    hour = date.strftime('%H')
    minutes = date.strftime('%M')
    
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

    return utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03

# Applying the Rayleigh correction
def applying_rayleigh_correction(utc_time, lons, lats, sun_zenith, data_ch01, data_ch02):
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
    # Adicionando os logos
    logo_noaa = plt.imread(dir_logos + 'NOAA_Logo.png')  # Lendo o arquivo do logo
    logo_goes = plt.imread(dir_logos + 'GOES_Logo.png')  # Lendo o arquivo do logo
    logo_cepagri = plt.imread(dir_logos + 'CEPAGRI-Logo.png')  # Lendo o arquivo do logo
    fig.figimage(logo_noaa, 32, 240, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_goes, 10, 160, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_cepagri, 10, 80, zorder=3, alpha=0.8, origin='upper')  # Plotando logo


def apply_cira_stretch(band_data):
    
    log_root = np.log10(0.0223)
    denom = (1.0 - log_root) * 0.75
    band_data *= 0.01
    band_data = band_data.clip(np.finfo(float).eps)
    band_data = np.log10(band_data)
    band_data -= log_root
    band_data /= denom
    return 1 + band_data


def delete_xml(dir_in):
    # Listar todos os arquivos no diretório
    arquivos = os.listdir(dir_in)

    # Filtrar os arquivos com extensão .nc.aux.xml
    arquivos_aux = [arquivo for arquivo in arquivos if arquivo.endswith('.nc.aux.xml')]

    # Excluir os arquivos auxiliares
    for arquivo in arquivos_aux:
        caminho_completo = os.path.join(dir_in, arquivo)
        os.remove(caminho_completo)


start = t.time()  
dir_in = f'/mnt/e/truecolor/'
#dir_in = f'/home/guimoura/Documentos/projeto/TrueColor/'
dir_shapefiles = f'{dir_in}shapefiles/'
dir_colortables = f'{dir_in}colortables/'
dir_logos = f'{dir_in}logos/'
dir_out = f'{dir_in}output/'

# Bands files
path_ch01 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C01_G16_s20232641220207_e20232641229515_c20232641229582.nc'
path_ch02 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C02_G16_s20232641220207_e20232641229515_c20232641229574.nc'
path_ch03 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C03_G16_s20232641220207_e20232641229515_c20232641229581.nc'

# path_ch01 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C01_G16_s20232632030204_e20232632039512_c20232632039570.nc'
# path_ch02 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C02_G16_s20232632030204_e20232632039512_c20232632039572.nc'
# path_ch03 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C03_G16_s20232632030204_e20232632039512_c20232632039562.nc'

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

# Resolution
resolution = 4.0

# Read the resolution
band_resolution_km = getattr(file_ch01, 'spatial_resolution')
band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

# Division factor to reduce image size
f = math.ceil(float(resolution / band_resolution_km))
# Limpando memória
del band_resolution_km

# Read the central longitude
longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin
# Read the semi major axis
a = file_ch01.variables['goes_imager_projection'].semi_major_axis
# Read the semi minor axis
b = file_ch01.variables['goes_imager_projection'].semi_minor_axis
# Calculate the image extent 
h = file_ch01.variables['goes_imager_projection'].perspective_point_height

# Pega data do file
add_seconds = int(file_ch01.variables['time_bounds'][0])
date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
date_formated = date.strftime('%Y-%m-%d %H:%M UTC')
date_file = date.strftime('%Y%m%d%H%M')

# Brasil
extent = [-90.0, -40.0, -20.0, 10.0]

# São Paulo
#extent = [-53.25, -26.0, -44.0, -19.5]  # Min lon, Min lat, Max lon, Max lat
# Choose the image resolution (the higher the number the faster the processing is)
#resolution = 1.0


#------------------------------------------------------------------------------------------------------
#-------------------------------------------Reprojetando----------------------------------------------#
#------------------------------------------------------------------------------------------------------
# band01
variable = "CMI"
# reprojetando
grid = remap(path_ch01, variable, extent, resolution, h, a, b, longitude)
# Lê o retorno da função
data_ch01 = grid.ReadAsArray()

#------------------------------------------------------------------------------------------------------
# band02
file_ch02 = Dataset(path_ch02)
variable = "CMI"
# reprojetando
grid = remap(path_ch02, variable, extent, resolution, h, a, b, longitude)
# Lê o retorno da função
data_ch02 = grid.ReadAsArray()

#------------------------------------------------------------------------------------------------------
# band03
file_ch03 = Dataset(path_ch03)
variable = "CMI"
# reprojetando
grid = remap(path_ch03, variable, extent, resolution, h, a, b, longitude)
# Lê o retorno da função 
data_ch03 = grid.ReadAsArray()
#------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------
# Calculating zenith correction
utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03 = calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03)

# Applying the Rayleigh correction
data_ch01, data_ch02 = applying_rayleigh_correction(utc_time, lons, lats, sun_zenith, data_ch01, data_ch02)

# Fazendo calculo True color
R = data_ch02
G = (data_ch01 + data_ch02) / 2 * 0.93 + 0.07 * data_ch03 
B = data_ch01

# Aplicando cira stretch
R = apply_cira_stretch(R)
G = apply_cira_stretch(G)
B = apply_cira_stretch(B)

# Create the RGB
RGB = np.stack([R, G, B], axis=2)		
#------------------------------------------------------------------------------------------------------

# Product Name
product = "truecolor" 

# Formatando a descricao a ser plotada na imagem
description = f' GOES-{satellite} Natural True Color {date_file}'
institution = "CEPAGRI - UNICAMP"

#------------------------------------------------------------------------------------------------------  
# Plot configuration
plot_config = {
"dpi": 150, 
"countries_color": 'cyan', "countries_width": data_ch01.shape[0] * 0.00100,
"continents_color": 'cyan', "continents_width": data_ch01.shape[0] * 0.00025,
"grid_color": 'white', "grid_width": data_ch01.shape[0] * 0.00025, "grid_interval": 10.0,
"title_text": "GOES-" + satellite[1:3] + " True color RGB ", "title_size": int(data_ch01.shape[1] * 0.005), "title_x_offset": int(data_ch01.shape[1] * 0.01), "title_y_offset": data_ch01.shape[0] - int(data_ch01.shape[0] * 0.016), 
"thick_interval": 0, "cbar_labelsize": int(data_ch01.shape[0] * 0.005), "cbar_labelpad": -int(data_ch01.shape[0] * 0.0),
"file_name_id_1": satellite,  "file_name_id_2": product
}

fig = plt.figure(figsize=(data_ch01.shape[1]/float(plot_config["dpi"] + 5), data_ch01.shape[0]/float(plot_config["dpi"])), dpi=plot_config["dpi"], frameon=True, edgecolor='black', facecolor='black')

# Define the projection
proj = ccrs.PlateCarree()

# Modificar caso a imagem esteja fora da área
ax = plt.axes([0, 0.02, 1, 1], projection=proj)
# Não sei o que muda :(
ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Adicionando o shapefile dos estados brasileiros
adicionando_shapefile('br', ax)

# Adicionando  linhas dos litorais
adicionando_linhas(ax)

# Define the image extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]

# Plot the image
img = ax.imshow(RGB, origin='upper', extent=img_extent, zorder=1)

# Adicionando descricao da imagem
adicionando_descricao_imagem(description, institution, ax, fig)

# Adicionando as logos
adicionando_logos(fig)

# Salvando a imagem de saida
plt.savefig(dir_out + plot_config["file_name_id_2"] + "_" + date_file + '.png', facecolor='black')

delete_xml(f'{dir_in}goes/')

print('Total processing time:', round((t.time() - start),2), 'seconds.') 