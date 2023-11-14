#!/opt/miniconda3/envs/goes/bin/python3
# -*- coding: utf-8 -*-
# Packages = conda create --name goes -c conda-forge matplotlib netcdf4 cartopy boto3 gdal scipy pandas scp
# Packages = apt install ffmpeg
# ======================================================================================================
# Manipulando imagens GOES-16 NetCDF's
# ===================================# Bibliotecas necessarias =========================================
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # Plotagem de resultados, textos, logos, etc.
import cartopy.crs as ccrs  # Utilitario para sistemas de referência e coordenadas
import cartopy.io.shapereader as shpreader  # Utilitario pa
from netCDF4 import Dataset  # Utilitario para a biblioteca NetCDF4
import numpy as np  # Suporte para arrays e matrizes multidimensionais, com diversas funções matemáticas para trabalhar com estas estruturas
import datetime  # Utilitario para datas e horas
import time  # Utilitario para trabalhar com tempos
import logging  # Utilitario para criar os logs
from shutil import copyfile  # Utilitario para copia de arquivos
from shapely.geometry import Point
from shapely.ops import cascaded_union
from dirs import get_dirs
from utilities import area_para_recorte
from utilities import adicionando_shapefile
from utilities import adicionando_linhas
from utilities import adicionando_descricao_imagem
from utilities import adicionando_logos
from utilities import apply_cira_stretch
from utilities import applying_rayleigh_correction
from utilities import calculating_lons_lats
from remap import remap
from datetime import timedelta, datetime
from utilities import download_prod
import os
import re
from multiprocessing import Process, Manager
from shapely.geometry import Point
import shapely.speedups
shapely.speedups.enable()

dirs = get_dirs()
dir_in = dirs['dir_in']
dir_out = dirs['dir_out']
dir_colortables = dirs['dir_colortables']
dir_maps = dirs['dir_maps']
dir_shapefiles = dirs['dir_shapefiles']
dir_temp = dirs['dir_temp']


def degrees(file_id):
        
        proj_info = file_id.variables['goes_imager_projection']
        lon_origin = proj_info.longitude_of_projection_origin
        H = proj_info.perspective_point_height + proj_info.semi_major_axis
        r_eq = proj_info.semi_major_axis
        r_pol = proj_info.semi_minor_axis

        # Data info
        lat_rad_1d = file_id.variables['x'][:]
        lon_rad_1d = file_id.variables['y'][:]

        # Create meshgrid filled with radian angles
        lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)

        # lat/lon calculus routine from satellite radian angle vectors
        lambda_0 = (lon_origin * np.pi) / 180.0

        a_var = np.power(np.sin(lat_rad), 2.0) + (np.power(np.cos(lat_rad), 2.0) * (
                np.power(np.cos(lon_rad), 2.0) + (((r_eq * r_eq) / (r_pol * r_pol)) * np.power(np.sin(lon_rad), 2.0))))
        b_var = -2.0 * H * np.cos(lat_rad) * np.cos(lon_rad)
        c_var = (H ** 2.0) - (r_eq ** 2.0)

        r_s = (-1.0 * b_var - np.sqrt((b_var ** 2) - (4.0 * a_var * c_var))) / (2.0 * a_var)

        s_x = r_s * np.cos(lat_rad) * np.cos(lon_rad)
        s_y = - r_s * np.sin(lat_rad)
        s_z = r_s * np.cos(lat_rad) * np.sin(lon_rad)

        lat = (180.0 / np.pi) * (np.arctan(((r_eq * r_eq) / (r_pol * r_pol)) * (s_z / np.sqrt(((H - s_x) * (H - s_x)) + (s_y * s_y)))))
        lon = (lambda_0 - np.arctan(s_y / (H - s_x))) * (180.0 / np.pi)

        return lat, lon


def save_txt(array, nome_arquivo_txt):
    # Checa se a matriz é vazia
    if len(array) == 0:
        print(f'{nome_arquivo_txt} vazia')
        pass
    else:
        # Criando nome do arquivo e diretório -- Mudar as barras para Linux
        with open(f"{dir_out}fdcf/{nome_arquivo_txt}.txt", 'w') as file:  # tomar cuidado se não ele fica criando diretorio
            for valor in array:
                valor = f"{valor[0]},{valor[1]}\n"
                file.write(valor)


def save_log_erro(array_errors, nome_arquivo_txt):
    # Checa se a matriz é vazia
    if len(array_errors) == 0:
        pass
    else:
        # Criando nome do arquivo e diretório -- Mudar as barras para Linux
        with open(f"{dir_out}fdcf/{nome_arquivo_txt}.txt", 'w') as file:
            for valor in array_errors:
                erro = f"{valor}\n"
                file.write(erro)
                                 

def download_arquivos_fdcf():
    
    start = time.time()
    
    # Desired extent
    product_name = "ABI-L2-FDCF"
    
    #Captura hora atual em UTC para download no site da Amazon
    data_hora_atual = datetime.utcnow()

    #Atrasa 10 min para entrar em conformidade com Amazon
    data_10_min = datetime.strftime(data_hora_atual-timedelta(minutes=10),'%Y%m%d%H%M')

    #Correção para poder fazer download em qualquer horário
    data_hora_download_file = data_10_min[0:11]+ '0'

    # 
    dir_fdcf = f'{dir_in}fdcf/'
    
    # Verifica se o dir existe
    os.makedirs(dir_fdcf, exist_ok=True)
    
    # Baixa o produto
    fdcf = download_prod(data_hora_download_file, product_name, dir_fdcf)
    
    #
    print(f'Tempo total {round((time.time() - start), 2)} segundos ')
    
    return fdcf


def fdcf_tabela_hot_spots(date, ax):    
    
    # Cria uma lista com os itens no diretorio temp que sao arquivos e se encaixa na expressao regular "^fdcf_.+_.+_br.txt$"
    fdcf_list = [name for name in os.listdir(f'{dir_out}fdcf') if os.path.isfile(os.path.join(f'{dir_out}fdcf', name)) and re.match(f'^fdcf_{date.strftime("%Y%m%d")}_.+_br.txt$', name)]

    # Cria matriz diária de pontos e log. de erros e faz o ‘loop’ nos arquivos do diretório.
    matriz_diaria = []

    ## Captura a lista com os "contornos" dos estados para separar o total de pontos individual
    geometry_estados = list(shpreader.Reader(dir_shapefiles + "divisao_estados/gadm40_BRA_1.shp").geometries())

    ## Os 26 estados estão em ordem alfabética no arquivo shapefile e o indice 27 armazena o total diário
    lista_totalpixels_uf = [0 for i in range(28)]
    log_erro = []
    for name in fdcf_list:
        try:
            with open(f'{dir_out}fdcf/{name}', 'r') as file:
                for linha in file.readlines():
                    linha = linha.split(',')
                    p = [float(linha[0]), float(linha[1])]  # Ponto em lat,lon
                    ## Checkagem para evitar contagem de pontos duplicados
                    if not (p in matriz_diaria):   
                        ax.plot(float(linha[1]), float(linha[0]), 'r+', ms=2.5, transform=ccrs.PlateCarree())  # plotando em lon,lat
                        matriz_diaria.append(p)
                        ## Contagem do total por estado
                        for j in range(27):
                            estado = geometry_estados[j]
                            if estado.covers(Point(float(linha[1]),float(linha[0]))):
                                lista_totalpixels_uf[j]+=1
                                break
                    else:
                        continue
        except Exception as erro:
            logging.info(f'Erro no processamento da matriz diária: {erro}')
            log_erro.append(f'{name} - Error: {erro}')
            pass

    # Calcula a soma do total de focos incendios diário
    lista_totalpixels_uf[27]+= np.sum(lista_totalpixels_uf)

    ### Manipulação dos dados para tabela
    lista_siglas_UF = ['AC','AL','AP','AM','BA','CE','DF','ES','GO','MA','MT','MS','MG','PA','PB','PR','PE','PI','RJ','RN',
                'RS','RO','RR','SC','SP','SE','TO','Total']
    columns = ('UF', 'Hot Spot')
    data = []
    for i in range(28):
        table_value = (lista_siglas_UF[i],lista_totalpixels_uf[i])
        data.append(table_value)
    
    ## Cria a tabela com os dados separados anteriormente
    tabela= matplotlib.table.table(ax =ax,cellText=data,colLabels=columns,cellLoc ='center',loc='lower right',rowLoc='center',
                            colLoc='center')
    tabela.auto_set_column_width(col=list(range(len(data))))

    # Cria os nomes dos arquivos diários e salva no "directory".

    logging.info("save_txt: matriz_diaria")
    save_txt(matriz_diaria, f'fdcf_{date.strftime("%Y%m%d")}_br')
    save_log_erro(log_erro, f'fdcf_{date.strftime("%Y%m%d")}_errors')

    return fdcf_list

             
def process_fdcf(fdcf, v_extent, fdcf_diario):
    
    global dir_colortables, dir_in, dir_out

    extent, resolution = area_para_recorte(v_extent)

    # Captura a hora para contagem do tempo de processamento da imagem
    processing_start_time = time.time()
    
    # Lendo imagem FDCF
    fire_mask = Dataset(fdcf)

    # Coletando do nome da imagem a data/hora
    dtime_fdcf = fdcf[fdcf.find("ABI-L2-FDCF-M6_G16_s") + 20:fdcf.find("_e") - 1]
    date = (datetime.strptime(dtime_fdcf, '%Y%j%H%M%S'))
    date_img = date.strftime('%d-%b-%Y')
    date_file = date.strftime('%Y%m%d_%H%M%S')
    
    matriz_pixels_fogo = []
    fire_mask_values = fire_mask.variables['Mask'][:, :]
    selected_fires = (fire_mask_values == 10) | (fire_mask_values == 11) | (fire_mask_values == 13) | (fire_mask_values == 30) | (fire_mask_values == 33)

    lat, lon = degrees(fire_mask)
    
    # separando Latitudes e Longitudes dos pontos
    p_lat = lat[selected_fires]
    p_lon = lon[selected_fires]
    
    brasil = list(shpreader.Reader(dir_shapefiles + "divisao_estados/gadm36_BRA_0.shp").geometries())
    
    def processar_parte(indice_inicio, indice_fim, p_lat_lista, p_lon, brasil, matriz_pixels_fogo):
        for i in range(indice_inicio, indice_fim):
            if i < len(p_lat_lista):
                if brasil[0].covers(Point(p_lon[i], p_lat_lista[i])):
                    p = (p_lat_lista[i], p_lon[i])
                    print(p)
                    matriz_pixels_fogo.append(p)


    num_partes = 12
    tamanho_parte = len(p_lat) // num_partes
    resto = len(p_lat) % num_partes

    # Dividir as partes com sobreposição
    partes = [range(i, i + tamanho_parte + 1) if i < resto else range(i, i + tamanho_parte) for i in range(0, len(p_lat), tamanho_parte)]

    # Criar uma lista compartilhada
    manager = Manager()
    matriz_pixels_fogo = manager.list()

    # Processar cada parte
    processos = []
    for parte_indices in partes:
        processo_parte = Process(target=processar_parte, args=(parte_indices.start, parte_indices.stop, p_lat.tolist(), p_lon, brasil, matriz_pixels_fogo))
        processos.append(processo_parte)
        processo_parte.start()

    # Aguardar a conclusão dos processos
    for processo in processos:
        processo.join()

    # Converter a lista compartilhada para uma lista padrão
    matriz_pixels_fogo = list(set(matriz_pixels_fogo))

    # Ordenar a matriz
    matriz_pixels_fogo.sort(reverse=True)
    print(matriz_pixels_fogo)
    save_txt(matriz_pixels_fogo, f'fdcf_{date.strftime("%Y%m%d_%H%M%S")}_br')

    # Le o arquivo de controle de quantidade de pontos
    try:
        with open(f'{dir_temp}band21_control.txt', 'r') as fo:
            control = fo.readline()
    except:
        control = 0

    # Verifica se as ocorrencias de pontos é maior que as anteriores, se sim, armazena a quantidade e as imagens para gerar fundo
    print("Len matriz_pixels_fogo: ", len(matriz_pixels_fogo), " int control: ", int(control))
    date_ini = datetime(date.year, date.month, date.day, int(13), int(00))
    date_end = datetime(date.year, date.month, date.day, int(18), int(1))
    
    if len(matriz_pixels_fogo) > int(control) and date_ini <= date <= date_end:
        with open(f'{dir_temp}band21_control.txt', 'w') as fo:
            fo.write(str(len(matriz_pixels_fogo)))

    if fdcf_diario:
        # Reiniciar contagem para verificar imagem com maior quantidade de pontos no dia
        with open(f'{dir_temp}band21_control.txt', 'w') as fo:
            fo.write(str(0))

        ch01 = f'{dir_in}fdcf/ch01.nc'
        ch02 = f'{dir_in}fdcf/ch02.nc'
        ch03 = f'{dir_in}fdcf/ch03.nc'
        
        file_ch01 = Dataset(ch01)

        # Lê a longitude central
        longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin

        # Lê a data do arquivo
        add_seconds = int(file_ch01.variables['time_bounds'][0])
        date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
        date_file = date.strftime('%Y%m%d_%H%M%S')
        date_img = date.strftime('%d-%b-%Y UTC')

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

        # Adicionando descricao da imagem.
        description = f"GOES-16 Natural True Color,    Fire Hot Spot em {date_img}"  # Esse espaço é necessário para adicionar o caractere na imagem
        institution = f'CEPAGRI - UNICAMP'

        # Definindo tamanho da imagem de saida.
        d_p_i = 150
        fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')
        
        # Utilizando projecao geoestacionaria no cartopy
        img_extent = [extent[0], extent[2], extent[1], extent[3]]
        ax = plt.axes(projection=ccrs.PlateCarree(), extent=img_extent)
        ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

        # Adicionando o shapefile dos estados brasileiros
        adicionando_shapefile(v_extent, ax)
        
        # Plotando a imagem RGB
        ax.imshow(RGB, origin='upper', extent=img_extent)

        # Adicionando descricao da imagem ##Apenas fire hot spot tem a cruz = true para descrição da imagem
        adicionando_descricao_imagem(description, institution, ax, fig, cruz=True)

        # Adicionando os logos
        adicionando_logos(fig)
        
        # Faz tabela de hot spots
        fdcf_list = fdcf_tabela_hot_spots(date, ax)
        
        for name in fdcf_list:
            os.remove(f'{dir_out}fdcf/{name}')
                
        # Salva imagem
        plt.savefig(f'{dir_out}fdcf/fdcf_{date_file}_br.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
        # Fecha a janela para limpar a memoria
        plt.close()

    # Realiza o log do calculo do tempo de processamento da imagem
    print(f'{fdcf} - {v_extent} - {str(round(time.time() - processing_start_time, 4))} segundos')

if __name__ == "__main__":
      
    fdcf = download_arquivos_fdcf()
    v_extent = 'br'
    fdcf_diario = False
    
    process_fdcf(f'{dir_in}fdcf/{fdcf}.nc', v_extent, fdcf_diario)
    os.remove(f'{dir_in}fdcf/{fdcf}.nc')