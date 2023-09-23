# Processamento de Imagens GOES-16

Este repositório contém código Python para processar imagens do satélite GOES-16 e criar imagens True Color e Air Mass. As imagens processadas são geradas a partir dos dados das bandas 01, 02 e 03 do GOES-16. As imagens True Color são imagens naturais da Terra, enquanto as imagens Air Mass representam a distribuição vertical de massas de ar na atmosfera.

## Estrutura de Pastas

O repositório está estruturado da seguinte forma:

```
/truecolor
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

Você também deve ter os dados de imagem das bandas 01, 02 e 03 do GOES-16 disponíveis no diretório de entrada especificado.

## Instalação

1. Clone este repositório para o seu ambiente local:

```bash
git clone https://github.com/seu-usuario/seu-repositorio.git
```

2. Instale as bibliotecas Python necessárias, se ainda não estiverem instaladas:

```bash
pip install matplotlib numpy netCDF4 cartopy pyorbital pyspectral
```

## Configuração

Antes de executar os scripts, você precisa configurar algumas variáveis no script `Processamento.py`:

- `dir_main`: O diretório raiz do projeto.
- `dir_in`: O diretório onde estão os dados de entrada das bandas 01, 02 e 03 do GOES-16.
- `dir_shapefiles`: O diretório onde estão os arquivos shapefile para adicionar informações geográficas.
- `dir_colortables`: O diretório onde estão as tabelas de cores.
- `dir_logos`: O diretório onde estão os logotipos.
- `dir_out`: O diretório onde as imagens processadas serão salvas.

Além disso, você pode configurar as variáveis `bands` e `process_br`/`process_sp` para controlar quais bandas e regiões você deseja processar.

## Uso

Para processar as imagens, execute o script `Processamento.py`. Ele chamará os scripts `truecolor.py` e `airmass_rgb.py` para processar as imagens True Color e Air Mass, respectivamente. Certifique-se de que seus dados de entrada estejam disponíveis no diretório especificado em `dir_in`.

```bash
python Processamento.py
```

As imagens processadas serão salvas no diretório especificado em `dir_out`.

## Sobre o Código

Os scripts `truecolor.py` e `airmass_rgb.py` realizam o processamento das imagens True Color e Air Mass, respectivamente. Eles incluem funções para reprojetar dados, aplicar correções atmosféricas e criar imagens finais.

## Licença

Este código é distribuído sob a licença MIT. Consulte o arquivo LICENSE para obter mais detalhes.

---