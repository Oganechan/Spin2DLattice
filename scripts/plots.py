import matplotlib.pyplot as plt
import numpy as np

def plot_simulation_data(filename):
    # Чтение данных
    data = np.loadtxt(filename)
    
    concentration = data[:, 0]
    temperature = data[:, 1]
    energy = data[:, 2]
    magnetization = data[:, 3]
    
    # Создание фигуры с двумя подграфиками
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # График энергии
    ax1.plot(temperature, energy, 'b-', marker='o', markersize=3, linewidth=1.5)
    ax1.set_ylabel('Энергия')
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f'Температурная зависимость энергии (c = {concentration[0]})')
    
    # График намагниченности
    ax2.plot(temperature, magnetization, 'r-', marker='s', markersize=3, linewidth=1.5)
    ax2.set_xlabel('Температура (K)')
    ax2.set_ylabel('Намагниченность')
    ax2.grid(True, alpha=0.3)
    ax2.set_title(f'Температурная зависимость намагниченности (c = {concentration[0]})')
    
    plt.tight_layout()
    plt.savefig('temperature_dependencies.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return temperature, energy, magnetization

# Использование
if __name__ == "__main__":
    t, e, m = plot_simulation_data('../data/output/statistics_L20_rectangular.dat')