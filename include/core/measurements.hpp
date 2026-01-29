#ifndef MEASUREMENTS_HPP
#define MEASUREMENTS_HPP

#include <vector>

struct Measurements {
    std::vector<double> energies_;
    std::vector<double> magnetizations_;
    std::vector<double> order_parameters_;

    std::vector<double> concentrations_;
    std::vector<double> temperatures_;

    std::vector<double> mean_energy_;
    std::vector<double> mean_energy_sq_;
    std::vector<double> mean_energy_4th_;

    std::vector<double> mean_magnetization_;
    std::vector<double> mean_magnetization_abs_;
    std::vector<double> mean_magnetization_sq_;
    std::vector<double> mean_magnetization_4th_;

    std::vector<double> mean_order_parameter_;
    std::vector<double> mean_order_parameter_sq_;
    std::vector<double> mean_order_parameter_4th_;

    std::vector<double> specific_heat_;
    std::vector<double> magnetic_susceptibility_;
    std::vector<double> order_susceptibility_;

    std::vector<double> energy_binder_;
    std::vector<double> magnetic_binder_;
    std::vector<double> order_binder_;
};

#endif // MEASUREMENTS_HPP