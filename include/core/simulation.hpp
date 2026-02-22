#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "../algorithm/swendsenwang.hpp"
#include "../lattice/atoms.hpp"
#include "config.hpp"
#include "data.hpp"
#include <memory>
#include <string>

class Simulation {
  public:
    explicit Simulation(const Config &config,
                        const std::string &base_output_dir);

    void run();

  private:
    lattice::Atoms atoms_;
    physics::Calculator calculator_;
    std::unique_ptr<algorithm::ISwendsenWang> swendsenwang_;
    Data data_;

    const std::string base_output_dir_;

    const int32_t number_measures_;
    const std::string scan_type_;
    const double scan_start_, scan_step_, scan_end_;

    void run_production();
    void run_single_simulation();
    void run_temperature_scan();
    void run_concentration_scan();

    void dynamic_thermolysis();
};

#endif // SIMULATION_HPP