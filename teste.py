from datetime import datetime, timedelta
from dirs import get_dirs
from utilities import download_prod

dirs = get_dirs()
# Importando dirs do modulo dirs.py

dir_in = dirs['dir_in']

#Captura hora atual em UTC para download no site da Amazon
data_hora_atual = datetime.utcnow()

#Atrasa 10 min para entrar em conformidade com Amazon
data_10_min = datetime.strftime(data_hora_atual-timedelta(hours=1),'%Y%m%d%H')

# Donwload do produto LST 2km Full disk
download_prod(data_10_min,'ABI-L2-LST2KMF',f'{dir_in}lst/')
