#include "simulation.hpp"
#include "config.hpp"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>

Simulation::Simulation(const Config &config,
                       const std::string &output_directory)
    : atoms_(config), calculator_(atoms_, config),
      metropolis_(atoms_, calculator_),
      output_directory_(
          output_directory + "/" + config.get<std::string>("material.name") +
          "/L" + std::to_string(config.get<int32_t>("lattice.system_size")) +
          "_N" +
          std::to_string(config.get<int32_t>("simulation.number_measures")) +
          ".dat"),
      data_(calculator_, output_directory_),
      number_measures_(config.get<int32_t>("simulation.number_measures")),
      scan_type_(config.get<std::string>("simulation.scan_type")),
      scan_start_(config.get<double>("simulation.scan_start")),
      scan_step_(config.get<double>("simulation.scan_step")),
      scan_end_(config.get<double>("simulation.scan_end")) {
    std::filesystem::create_directories(
        output_directory + "/" + config.get<std::string>("material.name"));
}

void Simulation::run() {
    auto start_time = std::chrono::steady_clock::now();

    std::ofstream file(output_directory_, std::ios::trunc);

    if (scan_type_ == "temperature")
        run_temperature_scan();
    else
        run_concentration_scan();

    auto end_time = std::chrono::steady_clock::now();

    std::cout << "Execution time: "
              << std::chrono::duration<double>(end_time - start_time).count()
              << std::endl;
}

void Simulation::run_single_simulation() {
    data_.reset();

    metropolis_.sweep(10 * atoms_.get_magnetic_count());

    for (int32_t measure = 0; measure < number_measures_; ++measure) {
        metropolis_.sweep();
        data_.measure();
    }

    data_.save();
}

void Simulation::run_temperature_scan() {
    atoms_.set_random_defects(1.0 - calculator_.get_concentration());
    atoms_.initialize_random();

    for (double T = scan_end_; T >= scan_start_ - 1e-5; T -= scan_step_) {
        calculator_.set_temperature(T);
        run_single_simulation();
    }
}

void Simulation::run_concentration_scan() {
    for (double c = scan_start_; c <= scan_end_ + 1e-5; c += scan_step_) {
        atoms_.set_random_defects(1.0 - c);
        run_single_simulation();
    }
}
