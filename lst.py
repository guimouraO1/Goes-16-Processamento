# Required modules
from netCDF4 import Dataset      
import matplotlib.pyplot as plt  # Plotting library
from dirs import get_dirs
from datetime import datetime, timedelta
from utilities import area_para_recorte, remap, adicionando_shapefile, adicionando_linhas, adicionando_descricao_imagem, adicionando_logos
import time
import cartopy, cartopy.crs as ccrs 
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
import logging


def process_lst(file, v_extent):

    # Captura a hora para contagem do tempo de processamento da imagem
    start = time.time()

    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)
    satellite = '16'
    variable = 'LST'
    
    # Reprojetando
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
    
    # Mask values less than -20 degrees
    data_lst = np.ma.masked_where((data_lst < min_temp), data_lst)

    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-{satellite} Land Surface Temperature (°C) {date_img}'
    institution = "CEPAGRI - UNICAMP"

    # Definindo tamanho da imagem de saida
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat
    
    # Criando imagem de fundo Natural Earth 1
    # raster = gdal.Open(f'{dir_maps}HYP_HR_SR_OB_DR.tif')
    # min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]
    # raster = gdal.Translate(f'{dir_maps}naturalEarth_sp.tif', raster, projWin = [min_lon, max_lat, max_lon, min_lat])
    
    if v_extent == 'sp':
        raster = gdal.Open(f'{dir_maps}naturalEarth_sp.tif')
    else:
        raster = gdal.Open(f'{dir_maps}naturalEarth_br.tif')
    
    # Lendo o RGB 
    array = raster.ReadAsArray()
    R = array[0,:,:].astype(float) / 255
    G = array[1,:,:].astype(float) / 255
    B = array[2,:,:].astype(float) / 255
    
    R[R==4] = 0
    G[G==5] = 0
    B[B==15] = 0
    rgb = np.stack([R, G, B], axis=2)
    
    # PLotando imagem de fundo
    ax.imshow(rgb, extent=img_extent)
        
    # Plotando a imagem Spectral_r
    img = ax.imshow(data_lst, origin='upper',vmin=-25, vmax=60, extent=img_extent, zorder=2, cmap='jet')

    # Adicionando barra da paleta de cores de acordo com o canal
    cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
    cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[-10, 10, 30, 50])
    cb.ax.set_xticklabels(['-10', '10', '30','50'])
    cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
    cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
    cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de coreses

    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)
    
    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)

    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)

    type_lst = 'lst'
    
    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}{type_lst}/{type_lst}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

    # Fecha a janela para limpar a memoria
    plt.close()
    
    print(f'Total processing time: {round((time.time() - start),2)} seconds.')


def iniciar_processo_lst(p_br, p_sp, bands, new_bands):
    # Checagem se e possivel gerar imagem Land Surface Temperature
    if bands['23']:
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS LAND SURFACE TEMPERATURE "BR"...')
            # Pega o nome do produto LST2KMF 
            lst = new_bands['23']
            # Pega o local do produto 
            file = f'{dir_in}lst/{lst}'
            try:
                # Inicia o Processamento
                process_lst(file, 'br')
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")
                
        if p_sp:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS LAND SURFACE TEMPERATURE "SP"...')
            try:
                # Inicia o Processamento
                process_lst(file, 'sp')
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")


if __name__ == "__main__":
    
    # Importa os diretórios do módulo dirs.py
    dirs = get_dirs()
    dir_in = dirs['dir_in']
    dir_out = dirs['dir_out']
    dir_colortables = dirs['dir_colortables']
    dir_maps = dirs['dir_maps']
    
    # Coloque as badas em goes/band0? e coloque o nome do arquivo aqui
    new_bands = {'23': 'OR_ABI-L2-LST2KMF-M6_G16_s20233000900210_e20233000909518_c20233000910381.nc'}

    # Cria um dicionário de dados airmass True para processar a imagem
    bands = {}
    bands['23'] = True

    # Define as Imagens que vão ser processadas
    p_br = True
    p_sp = True

    # Inicia a função de processamento 
    iniciar_processo_lst(p_br, p_sp, bands, new_bands)