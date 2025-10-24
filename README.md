# Spin2DLattice

A high-performance C++ implementation of Monte Carlo simulation for magnetic spin systems.

## Authors

- **Leonid Brykin** ([@d3adand3nd](https://github.com/d3adand3nd)) - Project inspiration and scientific guidance  
- **Murlyka** ([@murlyka-coder](https://github.com/murlyka-coder)) - Implementation and development

## Features

- **Heisenberg spin model** with 3D vector spins
- **Multiple crystal lattices**: rectangular, c_rectangular, triangular, honeycomb, kagome, lieb
- **Metropolis algorithm** for thermal equilibrium
- **Configurable parameters** via JSON

## Quick Start

### Clone and build
```bash
git clone https://github.com/Oganechan/Spin2DLattice
mkdir ./Spin2DLattice/build && cd ./Spin2DLattice/build
g++ ../include/core/data.cpp ../include/core/simulation.cpp ../include/lattice/atoms.cpp ../include/lattice/geometry.cpp ../include/physics/calculator.cpp ../include/algorithm/metropolis.cpp ../include/main.cpp -o Spim2DLattice
```

### Run simulation
```bash
./Lattice2DLattice
```

## Configuration example

### JSON Format

```json
{
    "lattice": {
        "crystal_type": "rectangular",
        "system_size": 20,
        "lattice_constant_a": 1.0,
        "lattice_constant_b": 1.0,
        "boundary_type": "periodic"
    },
    "physical": {
        "external_magnetic_field": [
            0.0,
            0.0,
            0.0
        ],
        "exchange_constants": [
            1.0,
            0.1,
            0.01
        ],
        "temperature": 25.0
    },
    "simulation": {
        "equilibration_sweeps": 1000,
        "production_sweeps": 10000,
        "measurement_interval": 10
    }
}
```
