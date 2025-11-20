#include "simulation.hpp"
#include "config.hpp"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>

Simulation::Simulation(const Config &config,
                       const std::string &output_directory)
    : output_directory_(output_directory),
      output_postfix_(
          "_L" + std::to_string(config.get<int32_t>("lattice.system_size")) +
          "_" + config.get<std::string>("lattice.crystal_type")),
      atoms_(config), calculator_(atoms_), metropolis_(atoms_),
      data_(output_directory),
      number_measures_(config.get<int32_t>("simulation.number_measures")) {}

void Simulation::run() {
    auto start_time = std::chrono::steady_clock::now();

    atoms_.initialize_random();
    run_production();

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "Execution time: "
              << std::chrono::duration<double>(end_time - start_time).count()
              << std::endl;
}

void Simulation::run_production() {
    data_.reset(output_postfix_);

    for (double c = 0.0; c <= 1.00; c += 0.01) {
        for (int measure = 0; measure < number_measures_; ++measure) {
            atoms_.set_random_defects(1.0 - c);
            for (int32_t sweep = 0; sweep < atoms_.get_magnetic_count();
                 ++sweep)
                metropolis_.sweep();

            data_.measure(calculator_);
        }

        data_.save(c);
    }
}
