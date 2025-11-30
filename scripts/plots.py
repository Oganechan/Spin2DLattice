import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath', type=str)
    parser.add_argument('--output', '-o',type=str, default='plot.png')
    
    args = parser.parse_args()

    data = np.loadtxt(args.filepath)

    concentration = data[:, 0]
    temperature = data[:, 1]
    energy = data[:, 2]
    magnetization = data[:, 3]

    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(8,10))

    ax1.plot(temperature,energy)
    ax1.set_xlabel('Temperature')
    ax1.set_ylabel('Energy')

    ax2.plot(temperature,magnetization)
    ax2.set_xlabel('Temperature')
    ax2.set_ylabel('Concentration')

    plt.tight_layout()
    plt.savefig(args.output)


if __name__ == "__main__":
    main()