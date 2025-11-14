#pragma once

#include "../algorithm/metropolis.h"
#include "../lattice/atoms.h"
#include "config.h"
#include "data.h"
#include <chrono>
#include <iostream>

class Simulation {
  public:
    Simulation(const Config &config, const std::string &output_directory);

    void run();

  private:
    lattice::Atoms atoms_;
    physics::Calculator calculator_;
    algorithm::Metropolis metropolis_;
    Data data_;

    std::string output_directory_;

    const int32_t equilibration_sweeps_;
    const int32_t production_sweeps_;
    const int32_t measurement_interval_;

    void initialize_system();
    void run_equilibration();
    void run_production();
};