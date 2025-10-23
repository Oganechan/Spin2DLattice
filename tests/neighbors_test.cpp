#include <iostream>
#include "tests.h"

void print_neighbors_indexes(const lattice::Geometry &geometry)
{
    const auto &neighbors = geometry.get_neighbor_table();
    const int32_t num_atoms = geometry.get_atom_count();
    const int32_t num_shells = geometry.get_shell_count();

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

    fill_test_config(config, 10);
    config.set("physical.exchange_constants", std::vector<double>{1.0, 0.25});

    lattice::Atoms atoms(config);
    lattice::Geometry geometry = atoms.get_geometry();
    print_neighbors_indexes(geometry);
}