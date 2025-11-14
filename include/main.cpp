#include "core/simulation.hpp"

int main() {
    Config config;
    config.load("../data/input/default.json");

    Simulation simulation(config, "../data/output");
    simulation.run();

    return 0;
}