import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filepath', type=str)
parser.add_argument('--output', '-o',type=str, default='lattice.png')
    
args = parser.parse_args()

df = pd.read_csv(args.filepath)

unique_sublattices = df['sublattice'].unique()
color_map = plt.cm.tab10 
sublattice_to_color = {sl: color_map(i/len(unique_sublattices)) 
                       for i, sl in enumerate(unique_sublattices)}

fig, ax = plt.subplots(figsize=(10, 8))
for sublattice, group in df.groupby('sublattice'):
    ax.scatter(group['x'], group['y'], 
               label=f'Sublattice {sublattice}',
               color=sublattice_to_color[sublattice],
               s=50, alpha=0.7, edgecolors='black', linewidth=0.5)

ax.set_xlabel('X coordinate', fontsize=12)
ax.set_ylabel('Y coordinate', fontsize=12)
ax.set_title(f'Atomic lattice visualization ({len(df)} atoms)', fontsize=14)
ax.grid(True, alpha=0.3, linestyle='--')
ax.legend(title='Sublattice')
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(args.output, dpi=300, bbox_inches='tight')