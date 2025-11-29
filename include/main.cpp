#include "core/simulation.hpp"
#include "utils/random.hpp"
#include <iostream>

int main() {
    Random::initialize_thread_based();

    try {
        Config config;
        config.load("../data/input/default.json");

        Simulation simulation(config, "../data/output");
        simulation.run();
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}