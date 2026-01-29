import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filepath', type=str)
parser.add_argument('--output', '-o',type=str, default='plot.png')
    
args = parser.parse_args()

data = pd.read_csv(args.filepath)

temperature = data['temperature']
energy = data['energy']
specific_heat = data['specific_heat']
energy_cumulant = data['energy_cumulant']
order_parameter = data['order_parameter']
order_susceptibility = data['susceptibility_P']
order_cumulant = data['order_cumulant']


fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6,1,figsize=(8,12))

ax1.plot(temperature, energy)
ax1.set_xlabel('Temperature')
ax1.set_ylabel('Energy')

ax2.plot(temperature,specific_heat)
ax2.set_xlabel('Temperature')
ax2.set_ylabel('Specific Heat')

ax3.plot(temperature,energy_cumulant)
ax3.set_xlabel('Temperature')
ax3.set_ylabel("Energy Cumulant")

ax4.plot(temperature,order_parameter)
ax4.set_xlabel('Temperature')
ax4.set_ylabel('Order Parameter')

ax5.plot(temperature, order_susceptibility)
ax5.set_xlabel('Temperature')
ax5.set_ylabel('Order Susceptibility')

ax6.plot(temperature, order_cumulant)
ax6.set_xlabel('Temperature')
ax6.set_ylabel('Order Cumulant')

plt.tight_layout()
plt.savefig(args.output)
