#include "data.h"

Data::Data(const std::string &output_path)
    : output_path_(output_path) {}

void Data::measure(const physics::Calculator &calculator)
{
    double energy = calculator.calculate_total_energy();
    double magnetization = calculator.calculate_total_magnetization();

    energies_.push_back(energy);
    magnetizations_.push_back(magnetization);
}

void Data::save()
{
    compute_statistics();
    save_time_series();
}

void Data::reset()
{
    energies_.clear();
    magnetizations_.clear();

    double mean_energy_ = 0.0;
    double mean_magnetization_ = 0.0;
}

void Data::compute_statistics()
{
    if (energies_.empty())
        return;

    double energy_sum = 0.0, magnetization_sum = 0.0;

    for (int32_t i = 0; i < energies_.size(); ++i)
    {
        energy_sum += energies_[i];
        magnetization_sum += magnetizations_[i];
    }

    mean_energy_ = energy_sum / energies_.size();
    mean_magnetization_ = magnetization_sum / magnetizations_.size();
}

void Data::save_time_series() const
{
    std::ofstream file(output_path_ + "/timeseries.dat");

    file << "#step energy magnetization\n";

    for (int32_t i = 0; i < energies_.size(); ++i)
        file << i << " " << energies_[i] << " " << magnetizations_[i] << "\n";
}