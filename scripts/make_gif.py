import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import argparse
import glob
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('dirpath', type=str)
parser.add_argument('-o', '--output', type=str, default='animation.gif')
parser.add_argument('--fps',type=int, default=30)
args = parser.parse_args()

snapshot_files = sorted(glob.glob(os.path.join(args.dirpath, 'snap_*.csv')))

fig, ax = plt.subplots(figsize=(8,8))

def extract_temperature(filename):
    match = re.search(r'T(\d+\.?\d*)', filename)
    return float(match.group(1))

def update(frame):
    ax.clear()
    df = pd.read_csv(snapshot_files[frame])
    colors = df['uz'].map({1: 'red', -1: 'blue'})
    scatter = ax.scatter(df['x'], df['y'], c=colors)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(df['x'].min() - 1.0, df['x'].max() + 1.0)
    ax.set_ylim(df['y'].min() - 1.0, df['y'].max() + 1.0)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    temp = extract_temperature(snapshot_files[frame])
    ax.set_title(f'T = {temp:.2f} K')
    
    return scatter

anim = animation.FuncAnimation(fig, update, frames=len(snapshot_files))
anim.save(args.output, fps=args.fps)

plt.close()