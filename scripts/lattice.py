import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib.collections import LineCollection

def draw_connections(ax, df, max_distance = 0.6):
    points = df[['x','y']].values
    connections = []

    for i in range(len(points)):
        for j in range(i+1, len(points)):
            dist = np.linalg.norm(points[i] - points[j])
            if dist < max_distance:
                connections.append([points[i], points[j]])
    
    if connections:
        lc = LineCollection(connections, colors='#404040', linewidths=1.2, alpha=0.35, zorder=1)
        ax.add_collection(lc)

def plot_lattice_plain(df, output_file):
    fig, ax = plt.subplots(figsize=(3.5, 6.5), dpi=10)
    draw_connections(ax,df)
    ax.scatter(df['x'], df['y'], color='#1a1a1a', s=45, alpha=1, zorder=2, edgecolors='none')
    ax.set_axis_off()
    ax.set_aspect('equal')
    margin = -0.7
    ax.set_xlim(df['x'].min() - margin, df['x'].max() + margin)
    ax.set_ylim(df['y'].min() - margin, df['y'].max() + margin)
    plt.tight_layout(pad=0.5)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.15)
    plt.close()

def plot_lattice_colored(df, output_file):
    fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=100)
    draw_connections(ax, df, max_distance=1.5)
    unique_sublattices = df['sublattice'].unique()
    if len(unique_sublattices) <= 12:
        color_map = plt.cm.tab10
    else:
        color_map = plt.cm.tab20
    
    for i, sublattice in enumerate(unique_sublattices):
        group = df[df['sublattice'] == sublattice]
        color = color_map(i / max(len(unique_sublattices) - 1, 1))
        ax.scatter(group['x'], group['y'], color=color, s=60, alpha=0.95, 
                   edgecolors='black', linewidth=0.8, zorder=2)
    
    ax.set_axis_off()
    ax.set_aspect('equal')
    margin = -0.7
    ax.set_xlim(df['x'].min() - margin, df['x'].max() + margin)
    ax.set_ylim(df['y'].min() - margin, df['y'].max() + margin)
    plt.tight_layout(pad=0.5)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.15)
    plt.close()

parser = argparse.ArgumentParser()
parser.add_argument('filepath', type=str)
parser.add_argument('--output', '-o', type=str, default='lattice')

args = parser.parse_args()

df = pd.read_csv(args.filepath)

base_name = args.output.rsplit('.', 1)[0] if '.' in args.output else args.output
plain_file = f"{base_name}_plain.png"
colored_file = f"{base_name}_colored.png"

plot_lattice_plain(df, plain_file)
plot_lattice_colored(df, colored_file)