import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Función para generar un degradado de color
def create_colormap():
    #colors = [(0.145, 0.941, 0.941), (0.149, 0.678, 0.89), (0.149, 0.424, 0.89), (0.675, 0.09, 0.851), (0.761, 0.02, 0.463)]  # RGB: Amarillo - Naranja - Rojo
    colors = [(0.0, 0.0, 1.0),    # Azul
          (0.0, 0.5, 0.0),    # Verde Oscuro
          (1.0, 1.0, 0.0),    # Amarillo
          (1.0, 0.647, 0.0),  # Naranja
          (1.0, 0.0, 0.0)]    # Rojo
    cmap_name = 'custom_colormap'
    return LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

def plot_color_map(x_values_reaction,y_values_energy,y_values_probability):
    # Crear el gráfico principal
    fig, ax1 = plt.subplots(figsize=(8, 6), dpi=300)

    # Crear el segundo eje y (derecho) con el degradado de color
    ax2 = ax1.twinx()
    cmap = create_colormap()
    image = np.column_stack([y_values_probability[::-1], y_values_probability[::-1]])
    extent = [x_values_reaction[0], x_values_reaction[-1], 0, 0.5]
    im = ax2.imshow(image, cmap=cmap, aspect='auto', extent=extent, alpha=0.5)

    # Configuraciones del segundo eje y (derecho)
    ax2.set_ylabel('') 
    ax2.set_yticks([])
    #ax2.set_ylabel(r"Tunneling Probability", fontsize=7)
    #ax2.tick_params('y', colors='black', labelsize=7)

    # Graficar la curva principal (energía) en el primer eje y (izquierdo)
    ax1.plot(x_values_reaction, y_values_energy, color='black')

    # Configuraciones del primer eje y (izquierdo)
    ax1.set_xlabel(r"$s$ (amu$^{1/2}$ Bohr)", fontsize=7)
    ax1.set_ylabel(r"$V_a^G$ (kcal mol$^{-1}$)", fontsize=7)
    ax1.tick_params('x', labelsize=7)  # Cambiar el tamaño de los valores del eje x
    ax1.tick_params('y', labelsize=7)  # Cambiar el tamaño de los valores del eje y izquierdo
    ax1.set_ylim(min(y_values_energy), max(y_values_energy))

    # Añadir barra de color (escala de colores) al lado derecho del gráfico
    cbar = fig.colorbar(im, ax=ax2, pad=0.1)
    cbar.ax.tick_params(labelsize=7)  # Cambiar el tamaño de los valores de la barra de color
    cbar.set_label('Tunneling Probability', fontsize=7)  # Añadir título con tamaño de letra 8 a la barra de colores

    # Mostrar el gráfico
    plt.subplots_adjust(bottom=0.15)
    plt.show()