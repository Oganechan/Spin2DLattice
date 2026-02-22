#include "simulation.hpp"
#include "config.hpp"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>

Simulation::Simulation(const Config &config, const std::string &base_output_dir)
    : atoms_(config), calculator_(atoms_, config),
      base_output_dir_(base_output_dir),
      data_(calculator_, config, base_output_dir),
      number_measures_(config.get<int32_t>("simulation.number_measures")),
      scan_type_(config.get<std::string>("simulation.scan_type")),
      scan_start_(config.get<double>("simulation.scan_start")),
      scan_step_(config.get<double>("simulation.scan_step")),
      scan_end_(config.get<double>("simulation.scan_end")) {
    uint32_t dir_count =
        (atoms_.get_spin_model() == lattice::SpinModel::ISING
             ? 1
             : (atoms_.get_spin_model() == lattice::SpinModel::XY ? 2 : 3));
    swendsenwang_ =
        algorithm::SwendsenWangFactory::create(dir_count, atoms_, calculator_);
}

void Simulation::run() {
    auto start_time = std::chrono::steady_clock::now();

    if (scan_type_ == "temperature")
        run_temperature_scan();
    else
        run_concentration_scan();

    data_.save_finale();
    auto end_time = std::chrono::steady_clock::now();

    std::cout << "Execution time: "
              << std::chrono::duration<double>(end_time - start_time).count()
              << std::endl;
}

void Simulation::run_single_simulation() {

    swendsenwang_->sweep(10);

    for (int32_t measure = 0; measure < number_measures_; ++measure) {
        swendsenwang_->sweep(10);
        data_.measure();
    }

    data_.save_statistics();
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
