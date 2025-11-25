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
      number_measures_(config.get<int32_t>("simulation.number_measures")),
      scan_type_(config.get<std::string>("simulation.scan_type")),
      fixed_temperature_(config.get<double>("physical.temperature")),
      fixed_concentration_(config.get<double>("physical.concentration")) {}

void Simulation::run() {
    auto start_time = std::chrono::steady_clock::now();

    atoms_.initialize_random();
    metropolis_.sweep();

    if (scan_type_ == "temperature")
        run_temperature_scan();
    else
        run_concentration_scan();

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "Execution time: "
              << std::chrono::duration<double>(end_time - start_time).count()
              << std::endl;
}

void Simulation::run_single_simulation(double concentration,
                                       double temperature) {
    data_.reset(output_postfix_);

    atoms_.set_random_defects(1 - concentration);
    metropolis_.set_temperature(temperature);

    for (int32_t measure = 0; measure < number_measures_; ++measure) {
        metropolis_.sweep();
        data_.measure(calculator_);
    }
    data_.save(concentration, temperature);
}

void Simulation::run_temperature_scan() {
    for (double T = 0.1; T <= 2.0; T += 0.1)
        run_single_simulation(fixed_concentration_, T);
}

void Simulation::run_concentration_scan() {
    for (double c = 0.0; c <= 1.0; c += 0.01)
        run_single_simulation(c, fixed_temperature_);
}
