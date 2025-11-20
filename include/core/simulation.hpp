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
    std::string output_postfix_;

    const int32_t number_measures_;

    void run_production();
};

#endif // SIMULATION_HPP