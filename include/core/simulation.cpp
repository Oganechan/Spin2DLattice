#include "simulation.hpp"
#include <chrono>
#include <iostream>

Simulation::Simulation(const Config &config,
                       const std::string &output_directory)
    : output_directory_(output_directory), atoms_(config), calculator_(atoms_),
      metropolis_(atoms_), data_(output_directory),
      equilibration_sweeps_(
          config.get<int32_t>("simulation.equilibration_sweeps")),
      production_sweeps_(config.get<int32_t>("simulation.production_sweeps")),
      measurement_interval_(
          config.get<int32_t>("simulation.measurement_interval")) {}

void Simulation::run() {
    auto start_time = std::chrono::steady_clock::now();

    initialize_system();
    run_equilibration();
    run_production();

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "Execution time: "
              << std::chrono::duration<double>(end_time - start_time).count()
              << std::endl;
}

void Simulation::initialize_system() {
    atoms_.initialize_random();
    atoms_.set_random_defects(0.5);
}

void Simulation::run_equilibration() {
    if (equilibration_sweeps_ <= 0)
        return;

    for (int32_t sweep = 0; sweep < equilibration_sweeps_; ++sweep) {
        metropolis_.sweep();
    }
}

void Simulation::run_production() {
    data_.reset();

    for (int32_t sweep = 0; sweep < production_sweeps_; ++sweep) {
        metropolis_.sweep();

        if (sweep % measurement_interval_ == 0)
            data_.measure(calculator_);
    }

    data_.save();
}
