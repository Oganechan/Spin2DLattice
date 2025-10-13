#include <iostream>
#include "../include/lattice/geometry.h"

void print_neighbors_indexes(const Lattice &lattice)
{
    const auto &neighbors = lattice.get_neighbors();
    const int32_t num_atoms = lattice.get_num_atoms();
    const int32_t num_shells = lattice.get_num_shells();

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
    config.set("lattice.crystal", "rectangular");
    config.set("lattice.boundary", "periodic");
    config.set("system.exchange_interaction", std::vector<double>{1.0, 0.25});

    Lattice lattice(config);
    print_neighbors_indexes(lattice);
}