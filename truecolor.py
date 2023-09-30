import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset                                  # Lê / Escrever arquivos NetCDF4
import datetime                                              # Biblioteca para trabalhar com datas
from datetime import timedelta                               # Biblioteca para converter dia juliano para dd-mm-yyyy
import matplotlib.pyplot as plt                              # Biblioteca para plotar gráficos
import matplotlib.colors                                     # Cores do Matplotlib
import numpy as np                                           # Computação científica com Python
import cartopy, cartopy.crs as ccrs                          # Plotar mapas
import time as t                                             # Acesso e conversão de tempo
from remap import remap                                      # Importar a função Remap
from multiprocessing import Process  # Utilitario para multiprocessamento
import logging
from utilities import area_para_recorte
from utilities import adicionando_shapefile
from utilities import adicionando_linhas
from utilities import adicionando_descricao_imagem
from utilities import adicionando_logos
from utilities import apply_cira_stretch
from utilities import applying_rayleigh_correction
from utilities import calculating_lons_lats
from dirs import get_dirs

def process_truecolor(rgb_type, v_extent, ch01=None, ch02=None, ch03=None):
    global dir_out
    # Calcula o tempo de processamento True Color
    start = t.time() 
    # Lê a imagem da banda 01
    file_ch01 = Dataset(ch01)
    # Lê o identificador do satélite
    satellite = '16'
    # Lê a longitude central
    longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin

    # Lê a data do arquivo
    add_seconds = int(file_ch01.variables['time_bounds'][0])
    date = datetime.datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
    date_file = date.strftime('%Y%m%d_%H%M%S')
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')

    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)
    variable = "CMI"
    
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
    
    # If zenith angle is greater than 85°, the composite pixel is zero
    RGB[sun_zenith > 85] = 0
    # Create the mask for the regions with zero
    mask = (RGB == [0.0,0.0,0.0]).all(axis=2)
    # Apply the mask to overwrite the pixels
    RGB[mask] = [0,0,0]
    #------------------------------------------------------------------------------------------------------

    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-{satellite} Natural True Color {date_img}'
    institution = "CEPAGRI - UNICAMP"

    # Definindo tamanho da imagem de saida
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)

    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)

    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat
   
    # Plotando a imagem
    ax.imshow(RGB, origin='upper', extent=img_extent,  zorder=1)
    
    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)
    
    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}{rgb_type}/{rgb_type}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
    
    # Fecha a janela para limpar a memoria
    plt.close()

    # Tempo de processamento True color
    logging.info('Total processing time True Color:', round((t.time() - start),2), 'seconds.') 


def iniciar_processo_truelocor(p_br, p_sp, bands, process_br, process_sp, new_bands):
    
    # Checagem se e possivel gerar imagem TrueColor
    if bands['17']:
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS TRUECOLOR "BR"...')
            # Pegando nome das bandas 01, 02, 03
            ch01 = new_bands['01']
            ch02 = new_bands['02']
            ch03 = new_bands['03']
            # Montando dicionario de argumentos
            kwargs = {'ch01': f'{dir_in}band01/{ch01}', 
                      'ch02': f'{dir_in}band02/{ch02}', 
                      'ch03': f'{dir_in}band03/{ch03}'}
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process = Process(target=process_truecolor, args=("truecolor", "br"), kwargs=kwargs)
                # Adiciona o processo na lista de controle do processamento paralelo
                process_br.append(process)
                # Inicia o processo
                process.start()
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_br:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa a lista de processos
        process_br.clear()
        
        # Se a variavel de controle de processamento sp for True, realiza o processamento
        if p_sp:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS TRUECOLOR "SP"...')
            # Pegando nome das bandas 01, 02, 03
            ch01 = new_bands['01']
            ch02 = new_bands['02']
            ch03 = new_bands['03']
            # Montando dicionario de argumentos
            kwargs = {'ch01': f'{dir_in}band01/{ch01}', 'ch02': f'{dir_in}band02/{ch02}', 
                      'ch03': f'{dir_in}band03/{ch03}'}
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process = Process(target=process_truecolor, args=("truecolor", "sp"), kwargs=kwargs)
                # Adiciona o processo na lista de controle do processamento paralelo
                process_sp.append(process)
                # Inicia o processo
                process.start()
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_sp:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa a lista de processos
        process_sp.clear()


dirs = get_dirs()
# Importando dirs do modulo dirs.py

dir_in = dirs['dir_in']
dir_out = dirs['dir_out']
dir_shapefiles = dirs['dir_shapefiles']
dir_colortables = dirs['dir_colortables']
dir_logos = dirs['dir_logos']
dir_out = dirs['dir_out']

bands = {}
bands['17'] = True
process_br = []
process_sp = []
p_br = True
p_sp = False

# Coloque as badas em goes/band0? e coloque o nome do arquivo aqui
new_bands = { '01': f'OR_ABI-L2-CMIPF-M6C01_G16_s20232711120209_e20232711129518_c20232711129578.nc', 
              '02': f'OR_ABI-L2-CMIPF-M6C02_G16_s20232711120209_e20232711129517_c20232711129568.nc',
              '03': f'OR_ABI-L2-CMIPF-M6C03_G16_s20232711120209_e20232711129517_c20232711129567.nc'}

iniciar_processo_truelocor(p_br, p_sp, bands, process_br, process_sp, new_bands)