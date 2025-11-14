#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "../algorithm/metropolis.hpp"
#include "../lattice/atoms.hpp"
#include "config.hpp"
#include "data.hpp"

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

#endif // SIMULATION_HPP