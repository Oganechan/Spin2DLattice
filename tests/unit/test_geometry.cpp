#include "../../include/lattice/geometry.hpp"
#include <cstdint>
#include <fstream>

int main() {
    Config config;
    config.load("../data/input/test.json");

    lattice::Geometry geometry(config);

    const auto &positions = geometry.get_atom_positions();
    const auto &sublattices = geometry.get_atom_sublattices();
    const auto &neighbor_table = geometry.get_neighbor_table();

    std::ofstream file1("lattice.csv");
    file1 << "atom_id,x,y,sublattice\n";
    for (int32_t i = 0; i < geometry.get_atom_count(); ++i) {
        file1 << i << "," << positions[i][0] << "," << positions[i][1] << ","
              << sublattices[i] << "\n";
    }
    file1.close();

    std::ofstream file2("neighbor.csv");
    file2 << "atom_id,shell,neighbor_id\n";
    for (int32_t i = 0; i < geometry.get_atom_count(); ++i)
        for (int32_t shell = 0; shell < geometry.get_shell_count(); ++shell)
            for (int32_t neighbor_id : neighbor_table[i][shell])
                file2 << i << "," << shell << "," << neighbor_id << "\n";
    file2.close();

    return 0;
}