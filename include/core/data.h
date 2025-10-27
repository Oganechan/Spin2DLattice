#pragma once

#include <string>
#include <fstream>
#include "../lattice/atoms.h"
#include "../physics/calculator.h"

class Data
{
public:
    Data(const std::string &output_path);

    void measure(const physics::Calculator &calculator);
    void save();
    void reset();

    inline double get_mean_energy() const { return mean_energy_; }
    inline double get_mean_magnetization() const { return mean_magnetization_; }

private:
    std::vector<double> energies_;
    std::vector<double> magnetizations_;
    std::vector<std::array<double, 3>> magnetization_vectors_;

    double mean_energy_ = 0.0;
    double mean_magnetization_ = 0.0;

    std::string output_path_;

    void compute_statistics();
    void save_time_series() const;
    void save_statistics() const;
};