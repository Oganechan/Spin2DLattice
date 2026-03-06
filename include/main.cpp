#include "../external/argparse/argparse.hpp"
#include "core/simulation.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("Spin_2D_Lattice");
    program.add_argument("config").default_value("../data/input/default.json");
    program.add_argument("-o", "--output").default_value("../data/output");

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return EXIT_FAILURE;
    }

    try {
        auto config_path = program.get<std::string>("config");
        auto output_path = program.get<std::string>("-o");

        Config config;
        config.load(config_path);

        Simulation simulation(config, output_path);
        simulation.run();
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}