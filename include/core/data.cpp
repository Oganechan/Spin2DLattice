#include "data.hpp"
#include <cstdint>
#include <fstream>

Data::Data(const physics::Calculator &calculator,
           const std::string &output_path)
    : calculator_(calculator), output_path_(output_path) {}

void Data::measure() {
    double energy = calculator_.calculate_total_energy();
    double magnetization = calculator_.calculate_total_magnetization();
    double order_parameter = calculator_.calculate_total_order();

    energies_.push_back(energy);
    magnetizations_.push_back(magnetization);
    order_parameters_.push_back(order_parameter);
}

void Data::save() {
    compute_statistics();

    std::ofstream file(output_path_, std::ios::app);

    file << calculator_.get_concentration() << "\t"
         << calculator_.get_temperature()
         << "\t"

         // EnergyI
         << mean_energy_ << "\t" << specific_heat_
         << "\t"

         // Order parameter
         << mean_order_parameter_ << "\t" << order_susceptibility_ << "\t"
         << order_binder_ << "\n";
}

void Data::reset() {
    energies_.clear();
    magnetizations_.clear();
    order_parameters_.clear();

    mean_energy_ = 0.0;
    mean_energy_sq = 0.0;

    mean_order_parameter_ = 0.0;
    mean_order_parameter_abs_ = 0.0;
    mean_order_parameter_sq_ = 0.0;
    mean_order_parameter_4th_ = 0.0;

    specific_heat_ = 0.0;
    order_susceptibility_ = 0.0;
    order_binder_ = 0.0;
}

void Data::compute_statistics() {
    if (energies_.empty() || magnetizations_.empty() ||
        order_parameters_.empty())
        return;

    int32_t N = energies_.size();
    int32_t magnetic_count = calculator_.get_atoms().get_magnetic_count();
    double temperature = calculator_.get_temperature();

    if (temperature < 1.0e-5)
        return;

    double energy_sum = 0.0, energy_sum_sq = 0.0;
    double order_sum = 0.0, order_abs_sum = 0.0, order_sq_sum = 0.0,
           order_4th_sum = 0.0;

    for (int32_t i = 0; i < N; ++i) {
        double e = energies_[i];
        energy_sum += e;
        energy_sum_sq += e * e;

        double phi = order_parameters_[i];
        double phi_abs = std::abs(phi);
        order_sum += phi;
        order_abs_sum += phi_abs;
        order_sq_sum += phi * phi;
        order_4th_sum += phi * phi * phi * phi;
    }

    mean_energy_ = energy_sum / N;
    mean_energy_sq = energy_sum_sq / N;

    mean_order_parameter_ = order_sum / N;
    mean_order_parameter_abs_ = order_abs_sum / N;
    mean_order_parameter_sq_ = order_sq_sum / N;
    mean_order_parameter_4th_ = order_4th_sum / N;

    specific_heat_ = (mean_energy_sq - mean_energy_ * mean_energy_) /
                     (magnetic_count * temperature * temperature);

    order_susceptibility_ = (mean_order_parameter_sq_ -
                             mean_order_parameter_ * mean_order_parameter_) /
                            (magnetic_count * temperature);

    order_binder_ =
        1.0 - mean_order_parameter_4th_ /
                  (3.0 * mean_order_parameter_sq_ * mean_order_parameter_sq_);
}