#include "data.hpp"
#include <fstream>

Data::Data(const physics::Calculator &calculator,
           const std::string &output_path)
    : calculator_(calculator), output_path_(output_path) {}

void Data::measure() {
    double energy = calculator_.calculate_total_energy();
    double magnetization = calculator_.calculate_total_magnetization();

    energies_.push_back(energy);
    magnetizations_.push_back(magnetization);
}

void Data::save() {
    compute_statistics();

    std::ofstream file(output_path_, std::ios::app);
    file << calculator_.get_concentration() << "\t"
         << calculator_.get_temperature() << "\t" << mean_energy_ << "\t"
         << mean_magnetization_ << "\n";
}

void Data::reset() {
    energies_.clear();
    magnetizations_.clear();

    mean_energy_ = 0.0;
    mean_magnetization_ = 0.0;
}

void Data::compute_statistics() {
    if (energies_.empty() || magnetizations_.empty())
        return;

    double energy_sum = 0.0, magnetization_sum = 0.0;

    for (int32_t i = 0; i < energies_.size(); ++i) {
        energy_sum += energies_[i];
        magnetization_sum += magnetizations_[i];
    }

    mean_energy_ = energy_sum / energies_.size();
    mean_magnetization_ = magnetization_sum / magnetizations_.size();
}