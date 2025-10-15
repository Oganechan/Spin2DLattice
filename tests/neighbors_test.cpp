#include <iostream>
#include "../include/lattice/atoms.h"

void print_neighbors_indexes(const lattice::Geometry &geometry)
{
    const auto &neighbors = geometry.get_neighbors();
    const int32_t num_atoms = geometry.get_num_atoms();
    const int32_t num_shells = geometry.get_num_shells();

    for (int32_t atom_idx = 0; atom_idx < num_atoms; ++atom_idx)
    {
        std::cout << "ATOM " << atom_idx << ":\n";
        const auto &atom_neighbors = neighbors[atom_idx];
        const int32_t shells_to_print = std::min(num_shells, static_cast<int32_t>(atom_neighbors.size()));

        for (int32_t shell = 0; shell < shells_to_print; ++shell)
        {
            const auto &shell_neighbors = atom_neighbors[shell];
            std::cout << "  Shell " << shell << " (" << shell_neighbors.size() << " neighbors): ";

            for (auto neighbor_idx : shell_neighbors)
                std::cout << neighbor_idx << " ";

            std::cout << "\n";
        }
        std::cout << "================\n";
    }
}

int main()
{
    Config config;

    config.set("lattice.linear_size", 10);
    config.set("lattice.norm_a", 1.0);
    config.set("lattice.norm_b", 1.0);
    config.set("lattice.crystal_type", "rectangular");
    config.set("lattice.boundary_conditions", "periodic");
    config.set("physical.exchange_interaction", std::vector<double>{1.0, 0.25});
    config.set("physical.spin_model", "ising");

    lattice::Atoms atoms(config);
    lattice::Geometry geometry = atoms.geometry();
    print_neighbors_indexes(geometry);
}