# Spin2DLattice

A high-performance C++ implementation of Monte Carlo simulation for magnetic spin systems.

## Authors

- **Leonid Afremov** (afremov.ll@dvfu.ru) - Scientific supervisor
- **Nikita Legkikh** ([@oganechan](https://github.com/Oganechan)) - Implementation and development

## Quick Start

### Clone and build
```bash
git clone https://github.com/Oganechan/Spin2DLattice
mkdir ./Spin2DLattice/build && cd ./Spin2DLattice/build
g++ ../include/core/data.cpp ../include/core/simulation.cpp ../include/lattice/atoms.cpp ../include/lattice/geometry.cpp ../include/physics/calculator.cpp ../include/main.cpp -o Spin2DLattice
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
        "system_size": 32,
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
            -1.0
        ],
        "anisotropy_constant": -1.0,
        "temperature": 1.0,
        "concentration": 1.0
    },
    "simulation": {
        "number_measures": 1000,
        "scan_type": "temperature",
        "scan_start": 0.01,
        "scan_step": 0.05,
        "scan_end": 1.0
    }
}
```
