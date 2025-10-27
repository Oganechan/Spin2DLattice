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
    save_statistics();
}

void Data::reset()
{
    energies_.clear();
    magnetizations_.clear();
    magnetization_vectors_.clear();

    mean_energy_ = 0.0;
    mean_magnetization_ = 0.0;
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
    file << "step\tenergy\tmagnetization\n";

    for (int32_t i = 0; i < energies_.size(); ++i)
        file << i << "\t" << energies_[i] << "\t" << magnetizations_[i] << "\n";
}

void Data::save_statistics() const
{
    std::ofstream file(output_path_ + "/statistics.dat");
    file << "quantity\tvalue\n"
         << "energy\t" << mean_energy_ << "\n"
         << "magnetization\t" << mean_magnetization_ << "\n ";
}