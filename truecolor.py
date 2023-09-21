import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset                                  # Lê / Escrever arquivos NetCDF4
import datetime                                              # Biblioteca para trabalhar com datas
from datetime import timedelta                               # Biblioteca para converter dia juliano para dd-mm-yyyy
import matplotlib.pyplot as plt                              # Biblioteca para plotar gráficos
import matplotlib.colors                                     # Cores do Matplotlib
import numpy as np                                           # Computação científica com Python
import cartopy, cartopy.crs as ccrs                          # Plotar mapas
import cartopy.io.shapereader as shpreader                   # Importar shapefiles
import time as t                                             # Acesso e conversão de tempo
import math                                                  # Importar math
import os                                                   # Interfaces do sistema operacional
from remap import remap                                      # Importar a função Remap
from pyorbital import astronomy
from pyspectral.rayleigh import Rayleigh     # Correção atmosférica no espectro visível 
from pyorbital.astronomy import get_alt_az
from pyorbital.orbital import get_observer_look
import logging

def adicionando_linhas(ax):
    # Adicionando linhas da costa
    ax.coastlines(resolution='10m', color='cyan', linewidth=0.5)
    # Adicionando linhas das fronteiras
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='cyan', linewidth=0.5)
    # Adicionando paralelos e meridianos
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
        # Adicionando o shapefile dos estados brasileiros e da cidade de Campinas
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='cyan', facecolor='none', linewidth=0.7)
        shapefile = list(shpreader.Reader(dir_shapefiles + 'campinas/campinas').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='yellow', facecolor='none', linewidth=1)

def adicionando_descricao_imagem(description, institution, ax, fig, cruz=False):
    # Criando novos eixos de acordo com a posição da imagem
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
    utc_time = datetime.datetime(int(year), int(month), int(day), int(hour), int(minutes))
    sun_zenith = np.zeros((data_ch01.shape[0], data_ch01.shape[1]))
    sun_zenith = astronomy.sun_zenith_angle(utc_time, lons, lats)

    # Aplicar a correção zenital do sol
    data_ch01 = (data_ch01)/(np.cos(np.deg2rad(sun_zenith)))
    data_ch02 = (data_ch02)/(np.cos(np.deg2rad(sun_zenith)))
    data_ch03 = (data_ch03)/(np.cos(np.deg2rad(sun_zenith)))

    return utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03

# Aplicar a correção de Rayleigh
def applying_rayleigh_correction(utc_time, lons, lats, sun_zenith, data_ch01, data_ch02):
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

def delete_xml(dir_in):
    # Listar todos os arquivos no diretório
    arquivos = os.listdir(dir_in)

    # Filtrar os arquivos com extensão .nc.aux.xml
    arquivos_aux = [arquivo for arquivo in arquivos if arquivo.endswith('.nc.aux.xml')]

    # Excluir os arquivos auxiliares
    for arquivo in arquivos_aux:
        caminho_completo = os.path.join(dir_in, arquivo)
        os.remove(caminho_completo)



dir_in = f'/mnt/e/truecolor/'
#dir_in = f'/home/guimoura/Documentos/projeto/TrueColor/'
dir_shapefiles = f'{dir_in}shapefiles/'
dir_colortables = f'{dir_in}colortables/'
dir_logos = f'{dir_in}logos/'
dir_out = f'{dir_in}output/'


start = t.time() 
# Arquivos das bandas
path_ch01 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C01_G16_s20232641220207_e20232641229515_c20232641229582.nc'
path_ch02 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C02_G16_s20232641220207_e20232641229515_c20232641229574.nc'
path_ch03 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C03_G16_s20232641220207_e20232641229515_c20232641229581.nc'

# path_ch01 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C01_G16_s20232632030204_e20232632039512_c20232632039570.nc'
# path_ch02 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C02_G16_s20232632030204_e20232632039512_c20232632039572.nc'
# path_ch03 = f'{dir_in}goes/OR_ABI-L2-CMIPF-M6C03_G16_s20232632030204_e20232632039512_c20232632039562.nc'


# Lê a imagem da banda 01
file_ch01 = Dataset(path_ch01)

# Lê o identificador do satélite
satellite = getattr(file_ch01, 'platform_ID')

# Lê o número da banda 
band = str(file_ch01.variables['band_id'][0]).zfill(2)

# Resolução
resolution = 3.0

# Lê a resolução da banda
band_resolution_km = getattr(file_ch01, 'spatial_resolution')
band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

# Fator de divisão para reduzir o tamanho da imagem
f = math.ceil(float(resolution / band_resolution_km))
# Limpando memória
del band_resolution_km

# Lê a longitude central
longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin
# Lê o semi-eixo maior
a = file_ch01.variables['goes_imager_projection'].semi_major_axis
# Lê o semi-eixo menor
b = file_ch01.variables['goes_imager_projection'].semi_minor_axis
# Calcular a extensão da imagem
h = file_ch01.variables['goes_imager_projection'].perspective_point_height

# Lê a data do arquivo
add_seconds = int(file_ch01.variables['time_bounds'][0])
date = datetime.datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
date_formated = date.strftime('%Y-%m-%d %H:%M UTC')
date_file = date.strftime('%Y%m%d_%H%M%S')
date_img = date.strftime('%d-%b-%Y %H:%M UTC')

# Extensão geográfica (Brasil)
extent = [-90.0, -40.0, -20.0, 10.0]

# Extensão geográfica (São Paulo)
#extent = [-53.25, -26.0, -44.0, -19.5]  # Min lon, Min lat, Max lon, Max lat
# Choose the image resolution (the higher the number the faster the processing is)
#resolution = 1.0


#------------------------------------------------------------------------------------------------------#
#-------------------------------------------Reprojetando----------------------------------------------#
#------------------------------------------------------------------------------------------------------#
# band01
variable = "CMI"
# reprojetando
grid = remap(path_ch01, variable, extent, resolution, h, a, b, longitude)
# Lê o retorno da função
data_ch01 = grid.ReadAsArray()

#------------------------------------------------------------------------------------------------------
# band02
file_ch02 = Dataset(path_ch02)
# reprojetando
grid = remap(path_ch02, variable, extent, resolution, h, a, b, longitude)
# Lê o retorno da função
data_ch02 = grid.ReadAsArray()

#------------------------------------------------------------------------------------------------------
# band03
file_ch03 = Dataset(path_ch03)
# reprojetando
grid = remap(path_ch03, variable, extent, resolution, h, a, b, longitude)
# Lê o retorno da função 
data_ch03 = grid.ReadAsArray()
#------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------
# Calculando correção zenith
utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03 = calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03)

# Aplicando a correção de Rayleigh
data_ch01, data_ch02 = applying_rayleigh_correction(utc_time, lons, lats, sun_zenith, data_ch01, data_ch02)

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

# Nome do produto
product = "truecolor" 

# Formatando a descricao a ser plotada na imagem
description = f' GOES-{satellite} Natural True Color {date_img}'
institution = "CEPAGRI - UNICAMP"

#------------------------------------------------------------------------------------------------------  
# Formato da imagem
dpi = 150
fig = plt.figure(figsize=(data_ch01.shape[1]/float(dpi + 5), data_ch01.shape[0]/float(dpi)), dpi=dpi, frameon=True, edgecolor='black', facecolor='black')

# Define a projeção
proj = ccrs.PlateCarree()

# Modificar caso a imagem esteja fora da área 
ax = plt.axes([0, 0.02, 1, 1], projection=proj)
# Não sei o que muda :(
ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Adicionando o shapefile dos estados brasileiros
adicionando_shapefile('br', ax)

# Adicionando  linhas dos litorais
adicionando_linhas(ax)

# Pega a extensão da imagem
img_extent = [extent[0], extent[2], extent[1], extent[3]]

# Plota a imagem
img = ax.imshow(RGB, origin='upper', extent=img_extent, zorder=1)

# Adicionando descricao da imagem
adicionando_descricao_imagem(description, institution, ax, fig)

# Adicionando as logos
adicionando_logos(fig)

# Salvando a imagem de saida
plt.savefig(dir_out + product + "_" + date_file + '.png', facecolor='black')

# Deletando arquivos .nc.aux.xml da lib gdal
delete_xml(f'{dir_in}goes/')

# Tempo de processamento True color
logging.info('Total processing time:', round((t.time() - start),2), 'seconds.') 