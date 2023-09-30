

# Modifique aqui os diretórios para uso
def get_dirs():
    #dir_main = f'/home/guimoura/Documentos/Goes-16-Processamento/'
    dir_main = f'/mnt/e/TrueColor/'
    dirs = {
        'dir_in': f'{dir_main}goes/',
        'dir_main': dir_main,
        'dir_out': dir_main + 'output/',
        'dir_shapefiles': dir_main + 'shapefiles/',
        'dir_colortables': dir_main + 'colortables/',
        'dir_logos': dir_main + 'logos/',
    }
    
    return dirs

