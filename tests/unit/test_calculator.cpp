#include "../../include/core/config.hpp"
#include "../../include/lattice/atoms.hpp"
#include "../../include/physics/calculator.hpp"
#include <iomanip>
#include <iostream>

void test() {
    Config config;
    config.load("../data/input/default.json");

    lattice::Atoms atoms(config);
    physics::Calculator calculator(atoms, config);

    for (int i = 0; i < 100; ++i) {
        std::cout << atoms.get_defect_count() << "\n";
        atoms.set_random_defects(1.0 - 0.5);
    }
}

int main() {
    Random::initialize_thread_based();
    std::cout << std::fixed << std::setprecision(6);

    try {
        test();

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}