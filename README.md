# Processamento de Imagens GOES-16

Este repositório contém código Python para processar imagens do satélite GOES-16 e criar imagens True Color e Air Mass. As imagens processadas são geradas a partir dos dados das bandas 01, 02 e 03 para True Color e 08, 10, 12 e 13 para Air Mass. As imagens True Color são imagens naturais da Terra, enquanto as imagens Air Mass representam a distribuição vertical de massas de ar na atmosfera.

## Estrutura de Pastas

O repositório está estruturado da seguinte forma:

```
/Goes-16-Processamento
    /colortables
    /logos
    /shapefiles
    airmass_rgb.py
    remap.py
    truecolor.py
    
```

- **colortables:** Contém tabelas de cores usadas no processamento.
- **logos:** Contém logotipos para adicionar às imagens geradas.
- **shapefiles:** Contém arquivos shapefile usados para adicionar informações geográficas às imagens processadas.
- **airmass_rgb.py:** Script para processamento de imagens Air Mass.
- **remap.py:** Script para reprojeção de dados.
- **truecolor.py:** Script para processamento de imagens True Color.
- **Processamento.py:** Script principal para iniciar o processamento.

## Pré-requisitos

Antes de executar os scripts, você precisa ter as seguintes bibliotecas Python instaladas:

- `matplotlib`
- `numpy`
- `netCDF4`
- `cartopy`
- `pyorbital`
- `pyspectral`
- `multiprocessing`
- `logging`

Você também deve ter os dados de imagem das bandas 01, 02 e 03 para True Color e 08, 10, 12 e 13 para Air Mass disponíveis no diretório de entrada especificado.

Ex: 
    goes/band08/OR_ABI-L2-CMIPF-M6C08_G16_s20232641020207_e20232641029515_c20232641029592.nc

## Instalação

1. Clone este repositório para o seu ambiente local:

```bash
git clone https://github.com/guimouraO1/TrueColor.git
```

2. Instale as bibliotecas Python necessárias, se ainda não estiverem instaladas:

```bash
pip install matplotlib numpy netCDF4 cartopy pyorbital pyspectral
```
Ou utilize o conda

```bash
conda create --name goes -c conda-forge matplotlib numpy netCDF4 cartopy pyorbital pyspectral
```

## Configuração

Antes de executar os scripts, você precisa criar as pastas goes/ e band?? que vc vai utilizar e os output/airmass e output/truecolor.

-Configure algumas variáveis no script:

- `dir_main`: O diretório raiz do projeto.
- `new_bands`: Coloque os arquivos netCDF4 das bandas dentro do dicionário de dados de banda correta.
Certifique-se de que seus dados de entrada estejam disponíveis no diretório especificado em `dir_in`.

Além disso, você pode configurar as variáveis `bands` e `process_br`/`process_sp` para controlar quais bandas e regiões você deseja processar.

## Uso

Para processar as imagens, execute o script desejado `truecolor.py` ou `airmass_rgb.py` para processar as imagens True Color ou Air Mass. 

Truecolor
```bash
python truecolor.py
```
Air Mass
```bash
python airmass_rgb.py
```
As imagens processadas serão salvas no diretório especificado em `dir_out`.

## Sobre o Código

Os scripts `truecolor.py` e `airmass_rgb.py` realizam o processamento das imagens True Color e Air Mass, respectivamente. Eles incluem funções para reprojetar dados, aplicar correções atmosféricas e criar imagens finais.

## Licença

Este código é distribuído sob a licença MIT. Consulte o arquivo LICENSE para obter mais detalhes.

---