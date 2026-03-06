# Spin2DLattice

A high-performance C++ implementation of Monte Carlo simulation for magnetic spin systems.

## Authors

- **Leonid Afremov** (afremov.ll@dvfu.ru) - Scientific supervisor
- **Nikita Legkikh** ([@oganechan](https://github.com/Oganechan)) - Implementation and development

## Quick Start

### Clone and build
```bash
git clone https://github.com/Oganechan/Spin2DLattice
mkdir Spin2DLattice/build && cd Spin2DLattice/build
cmake -S .. -B . && make
```

### Run simulation
```bash
./Spin2DLattice <config_path> -o <output_directory_path>
```

### Plots builing
```bash
bash ./script/build_plots.sh <data_directory_path>
```

## Configuration example

### JSON Format

```json
{
    "lattice": {
        "system_size": 16,
        "lattice_constant_a": 1.0,
        "lattice_constant_b": 1.0,
        "crystal_type": "rectangular",
        "boundary_type": "periodic"
    },
    "physical": {
        "model_type": "Ising",
        "external_magnetic_field": [
            0.0,
            0.0,
            0.0
        ],
        "exchange_constants": [
            1.0
        ],
        "anisotropy_constant": 0.0,
        "temperature": 1.0,
        "concentration": 1.0
    },
    "simulation": {
        "number_measures": 10000,
        "number_therms": 50000
        "scan_type": "temperature",
        "scan_start": 0.0,
        "scan_step": 0.01,
        "scan_end": 5.0
    }
}
```
