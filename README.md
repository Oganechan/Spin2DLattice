# Spin2DLattice

A high-performance C++ implementation of Monte Carlo simulation for magnetic spin systems.

## Authors

- **Leonid Brykin** ([@d3adand3nd](https://github.com/d3adand3nd)) - Project inspiration and scientific guidance
- **Nikita Legkikh** ([@oganechan](https://github.com/Oganechan)) - Implementation and development
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
g++ ../include/core/data.cpp ../include/core/simulation.cpp ../include/lattice/atoms.cpp ../include/lattice/geometry.cpp ../include/physics/calculator.cpp ../include/algorithm/metropolis.cpp ../include/main.cpp -o Spin2DLattice
```

### Run simulation
```bash
./Spin2DLattice
```

## Configuration example

### JSON Format

```json
{
    "lattice": {
        "crystal_type": "triangular",
        "system_size": 10,
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
            1.0
        ],
        "temperature": 0.5
    },
    "simulation": {
        "number_measures": 100
    }
}
```
