#------------------------------------------------------------------------------------------------------
# Required modules
import matplotlib.pyplot as plt              # Import the Matplotlib package
import numpy as np                           # Import the Numpy package
from netCDF4 import Dataset                  # Import the NetCDF Python interface
import cartopy, cartopy.crs as ccrs          # Plot maps
from datetime import datetime, timedelta     # Library to convert julian day to dd-mm-yyyy
from truecolor import area_para_recorte
from truecolor import adicionando_shapefile
from truecolor import adicionando_linhas
from truecolor import adicionando_descricao_imagem
from truecolor import adicionando_logos
from remap import remap
import time as time    
import logging
from multiprocessing import Process  # Utilitario para multiprocessamento
from matplotlib.colors import LinearSegmentedColormap, to_rgba

###########################################################################
#              Script de Processamento para Air Mass Goes-16              #
###########################################################################
#  Metodo: De 10 em 10 minutos                                            #
#  Descricao: Processa imagens netCDF4 para criar imagem Air Mass         #
#  Autor: Guilherme de Moura Oliveira  <guimoura@unicamp.br>              #
#  Data: 22/09/2023                                                       #
#  Atualizacao: 22/09/2023                                                #
###########################################################################


def process_airmass(rgb_type, v_extent, path_ch08=None, path_ch10=None, path_ch12=None, path_ch13=None):
    global dir_out
    start = time.time()  

    # Read the file using the NetCDF library
    file_ch08 = Dataset(path_ch08)

    # Lê o identificador do satélite
    satellite = getattr(file_ch08, 'platform_ID')
    variable = "CMI"
    
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)

    # Read the resolution
    band_resolution_km = getattr(file_ch08, 'spatial_resolution')
    band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

    # Getting the file time and date
    add_seconds = int(file_ch08.variables['time_bounds'][0])
    date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
    date_file = date.strftime('%Y%m%d_%H%M%S')
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')


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

    # Definindo tamanho da imagem de saida
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)

    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)    
    
    # Defina as cores da colorbar
    colors = ['#b62007','#6f008b', '#0a0a8e', '#538234','#5C8C3A','#335a25', '#704c02', '#b57350', '#ffffff']

    # Crie uma lista de posições normalizadas para as cores
    color_positions = np.linspace(0, 1, len(colors))

    # Crie um dicionário de cores segmentadas
    cmap_dict = {'red': [], 'green': [], 'blue': [], 'alpha': []}

    for color in colors:
        rgba = to_rgba(color)  # Use a função to_rgba para converter a cor
        cmap_dict['red'].append((color_positions[colors.index(color)], rgba[0], rgba[0]))
        cmap_dict['green'].append((color_positions[colors.index(color)], rgba[1], rgba[1]))
        cmap_dict['blue'].append((color_positions[colors.index(color)], rgba[2], rgba[2]))
        cmap_dict['alpha'].append((color_positions[colors.index(color)], rgba[3], rgba[3]))

    # Crie a paleta de cores personalizada
    custom_cmap = LinearSegmentedColormap('CustomCmap', cmap_dict)
    
    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]

    # Plot the image
    img = ax.imshow(RGB, origin='upper', cmap=custom_cmap, extent=img_extent)
        
    # Adicionando barra da paleta de cores de acordo com o canal
    cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
    cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[0.2, 0.4, 0.6, 0.8])
    cb.ax.set_xticklabels(['0.2', '0.4', '0.6','0.8'])
    cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
    cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
    cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de coreses
    
    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)

    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}{rgb_type}/{rgb_type}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

    # Fecha a janela para limpar a memoria
    plt.close()

    # Tempo de processamento True color
    logging.info('Total processing time Airmass:', round((time.time() - start),2), 'seconds.') 


def iniciar_processo_truelocor(p_br, p_sp, bands, process_br, process_sp, new_bands):
    # Checagem se e possivel gerar imagem Air Mass
    if bands['22']:
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS AIR MASS "BR"...')
            # Pegando nome das bandas 08, 10, 12, 13
            ch08 = new_bands['08']
            ch10 = new_bands['10']
            ch12 = new_bands['12']
            ch13 = new_bands['13']
            
            # Montando dicionario de argumentos
            kwargs = {'path_ch08': f'{dir_in}band08/{ch08}', 
                      'path_ch10': f'{dir_in}band10/{ch10}', 
                      'path_ch12': f'{dir_in}band12/{ch12}',
                      'path_ch13': f'{dir_in}band13/{ch13}'
                      }
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process = Process(target=process_airmass, args=("airmass", "br"), kwargs=kwargs)
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
                        # Pegando nome das bandas 08, 10, 12, 13
            ch08 = new_bands['08']
            ch10 = new_bands['10']
            ch12 = new_bands['12']
            ch13 = new_bands['13']
            
            # Montando dicionario de argumentos
            kwargs = {'path_ch08': f'{dir_in}band08/{ch08}', 
                      'path_ch10': f'{dir_in}band10/{ch10}', 
                      'path_ch12': f'{dir_in}band12/{ch12}',
                      'path_ch13': f'{dir_in}band13/{ch13}'
                      }
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process = Process(target=process_airmass, args=("airmass", "sp"), kwargs=kwargs)
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



#dir_main = f'/home/guimoura/Documentos/Goes-16-Processamento/'
dir_main =  f'/mnt/e/truecolor/' 
dir_out = f'{dir_main}output/'
dir_in = f'{dir_main}goes/'
dir_shapefiles = f'{dir_main}shapefiles/'
dir_colortables = f'{dir_main}colortables/'
dir_logos = f'{dir_main}logos/'

bands = {}
bands['22'] = True
process_br = []
process_sp = []
p_br = True
p_sp = True

# Coloque as badas em goes/band0? e coloque o nome do arquivo aqui
new_bands = { '08': f'OR_ABI-L2-CMIPF-M6C08_G16_s20232641020207_e20232641029515_c20232641029592.nc', 
              '10': f'OR_ABI-L2-CMIPF-M6C10_G16_s20232641020207_e20232641029526_c20232641029583.nc',
              '12': f'OR_ABI-L2-CMIPF-M6C12_G16_s20232641020207_e20232641029521_c20232641029598.nc',
              '13': f'OR_ABI-L2-CMIPF-M6C13_G16_s20232641020207_e20232641029526_c20232641029588.nc'
              }

iniciar_processo_truelocor(p_br, p_sp, bands, process_br, process_sp, new_bands)