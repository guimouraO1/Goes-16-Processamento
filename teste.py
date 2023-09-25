import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba

# Defina as cores da colorbar
colors = ['#6f008b', '#b62007', '#0a0a8e', '#538234', '#335a25', '#704c02', '#b57350', '#ffffff']

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

# Crie uma imagem de exemplo usando a paleta de cores personalizada
data = np.random.rand(10, 10)  # Substitua isso pelos seus próprios dados
plt.imshow(data, cmap=custom_cmap)

# Crie a colorbar
cbar = plt.colorbar(orientation='horizontal', label='Valor', ticks=color_positions)
cbar.set_ticklabels([])  # Remova os rótulos dos marcadores da colorbar

# Exiba o gráfico
plt.show()
