#include "data.hpp"
#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>

Data::Data(const physics::Calculator &calculator, const Config &config,
           const std::string &base_output_dir)
    : calculator_(calculator), config_(config) {
    simulation_name_ = generate_simulation_name();
    output_dir_ = fs::path(base_output_dir) / simulation_name_;

    initialize();
}

void Data::measure() {
    double energy = calculator_.calculate_total_energy();
    double magnetization = calculator_.calculate_total_magnetization();
    double order_parameter = calculator_.calculate_total_order();

    measurements_.energies_.push_back(energy);
    measurements_.magnetizations_.push_back(magnetization);
    measurements_.order_parameters_.push_back(order_parameter);
}

void Data::reset() {
    measurements_.energies_.clear();
    measurements_.magnetizations_.clear();
    measurements_.order_parameters_.clear();
}

void Data::compute_statistics() {
    if (measurements_.energies_.empty())
        return;

    int32_t measurement_count = measurements_.energies_.size();
    int32_t magnetic_count = calculator_.get_atoms().get_magnetic_count();

    double T = calculator_.get_temperature();
    int32_t N = calculator_.get_atoms().get_magnetic_count();

    if (T == 0.0)
        return;

    double mean_energy = 0.0;
    double mean_energy_sq = 0.0;
    double mean_energy_4th = 0.0;

    double mean_magnetization = 0.0;
    double mean_magnetization_abs = 0.0;
    double mean_magnetization_sq = 0.0;
    double mean_magnetization_4th = 0.0;

    double mean_order_parameter = 0.0;
    double mean_order_parameter_sq = 0.0;
    double mean_order_parameter_4th = 0.0;

    for (int32_t i = 0; i < measurement_count; ++i) {
        double e = measurements_.energies_[i];
        mean_energy += e;
        mean_energy_sq += e * e;
        mean_energy_4th += e * e * e * e;

        double m = measurements_.magnetizations_[i];
        mean_magnetization += m;
        mean_magnetization_abs += std::abs(m);
        mean_magnetization_sq += m * m;
        mean_magnetization_4th += m * m * m * m;

        double p = measurements_.order_parameters_[i];
        mean_order_parameter += p;
        mean_order_parameter_sq += p * p;
        mean_order_parameter_4th += p * p * p * p;
    }

    mean_energy /= measurement_count;
    mean_energy_sq /= measurement_count;
    mean_energy_4th /= measurement_count;
    mean_magnetization /= measurement_count;
    mean_magnetization_abs /= measurement_count;
    mean_magnetization_sq /= measurement_count;
    mean_magnetization_4th /= measurement_count;
    mean_order_parameter /= measurement_count;
    mean_order_parameter_sq /= measurement_count;
    mean_order_parameter_4th /= measurement_count;

    double specific_heat =
        N * (mean_energy_sq - mean_energy * mean_energy) / (T * T);
    double energy_binder =
        1.0 - mean_energy_4th / (3.0 * mean_energy_sq * mean_energy_sq);

    double magnetic_susceptibility =
        N * (mean_magnetization_sq - mean_magnetization * mean_magnetization) /
        T;
    double magnetic_binder =
        1.0 - mean_magnetization_4th /
                  (3.0 * mean_magnetization_sq * mean_magnetization_sq);

    double order_susceptibility =
        N *
        (mean_order_parameter_sq -
         mean_order_parameter * mean_order_parameter) /
        T;
    double order_binder =
        1.0 - mean_order_parameter_4th /
                  (3.0 * mean_order_parameter_sq * mean_order_parameter_sq);

    measurements_.concentrations_.push_back(calculator_.get_concentration());
    measurements_.temperatures_.push_back(T);
    measurements_.mean_energy_.push_back(mean_energy);
    measurements_.mean_energy_sq_.push_back(mean_energy_sq);
    measurements_.mean_energy_4th_.push_back(mean_energy_4th);
    measurements_.mean_magnetization_.push_back(mean_magnetization);
    measurements_.mean_magnetization_abs_.push_back(mean_magnetization_abs);
    measurements_.mean_magnetization_sq_.push_back(mean_magnetization_sq);
    measurements_.mean_magnetization_4th_.push_back(mean_magnetization_4th);
    measurements_.mean_order_parameter_.push_back(mean_order_parameter);
    measurements_.mean_order_parameter_sq_.push_back(mean_order_parameter_sq);
    measurements_.mean_order_parameter_4th_.push_back(mean_order_parameter_4th);
    measurements_.specific_heat_.push_back(specific_heat);
    measurements_.energy_binder_.push_back(energy_binder);
    measurements_.magnetic_susceptibility_.push_back(magnetic_susceptibility);
    measurements_.magnetic_binder_.push_back(magnetic_binder);
    measurements_.order_susceptibility_.push_back(order_susceptibility);
    measurements_.order_binder_.push_back(order_binder);
}

void Data::save_statistics() {
    compute_statistics();

    std::ofstream file(
        output_dir_ / "measurements" /
        ("T_" + std::to_string(calculator_.get_temperature()) + ".csv"));

    file << "energy,magnetization,order_parameter\n";
    for (int i = 0; i < measurements_.energies_.size(); ++i)
        file << measurements_.energies_[i] << ","
             << measurements_.magnetizations_[i] << ","
             << measurements_.order_parameters_[i] << "\n";
    file.close();

    reset();
}

void Data::save_finale() {
    std::ofstream file(output_dir_ / "results.csv");

    file << "concentration,temperature,energy,magnetization_abs,order_"
            "parameter,specific_heat,susceptibility_M,susceptibility_P,energy_"
            "cumulant,magnetic_cumulant,order_cumulant\n";
    for (int32_t i = 0; i < measurements_.temperatures_.size(); ++i)
        file << format_double(measurements_.concentrations_[i]) << ","
             << format_double(measurements_.temperatures_[i]) << ","
             << format_double(measurements_.mean_energy_[i]) << ","
             << format_double(measurements_.mean_magnetization_abs_[i]) << ","
             << format_double(measurements_.mean_order_parameter_[i]) << ","
             << format_double(measurements_.specific_heat_[i]) << ","
             << format_double(measurements_.magnetic_susceptibility_[i]) << ","
             << format_double(measurements_.order_susceptibility_[i]) << ","
             << format_double(measurements_.energy_binder_[i]) << ","
             << format_double(measurements_.magnetic_binder_[i]) << ","
             << format_double(measurements_.order_binder_[i]) << "\n";
    file.close();
}

std::string Data::generate_simulation_name() const {
    std::ostringstream name;

    std::string crystal_type = config_.get<std::string>("lattice.crystal_type");
    name << crystal_type << "_";

    std::string size =
        std::to_string(config_.get<int32_t>("lattice.system_size"));
    name << size << "x" << size << "_";

    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm *timeinfo = std::localtime(&t);
    std::ostringstream timestamp;
    timestamp << std::put_time(timeinfo, "%Y%m%d_%H%M%S");

    name << timestamp.str();

    return name.str();
}

void Data::initialize() {
    fs::create_directories(output_dir_);
    create_directory_structure();

    save_config();
    save_geometry();
    save_neighbor_table();
}

void Data::create_directory_structure() {
    fs::create_directories(output_dir_ / "geometry");
    fs::create_directories(output_dir_ / "snapshots");
    fs::create_directories(output_dir_ / "measurements");
    fs::create_directories(output_dir_ / "logs");
}

void Data::save_config() const {
    fs::path config_file = output_dir_ / "config_backup.json";
    std::ofstream file(config_file);

    file << config_.get().dump(4);
    file.close();
}

void Data::save_geometry() const {
    fs::path atoms_file = output_dir_ / "geometry" / "atoms.csv";
    std::ofstream file(atoms_file);
    file << "atom_id,x,y,sublattice\n";

    const auto &geometry = calculator_.get_atoms().get_geometry();
    const auto &atoms = calculator_.get_atoms();

    const auto &positions = geometry.get_atom_positions();
    const auto &sublattices = geometry.get_atom_sublattices();

    for (int32_t i = 0; i < geometry.get_atom_count(); ++i) {
        file << i << "," << positions[i][0] << "," << positions[i][1] << ","
             << sublattices[i] << "\n";
    }
    file.close();
}

void Data::save_neighbor_table() const {
    fs::path neighbors_file = output_dir_ / "geometry" / "neighbors.csv";
    std::ofstream file(neighbors_file);
    file << "atom_id,shell,neighbor_id\n";

    const auto &geometry = calculator_.get_atoms().get_geometry();
    const auto &neighbor_table = geometry.get_neighbor_table();

    for (int32_t i = 0; i < geometry.get_atom_count(); ++i)
        for (int32_t shell = 0; shell < geometry.get_shell_count(); ++shell)
            for (int32_t neighbor_id : neighbor_table[i][shell])
                file << i << "," << shell << "," << neighbor_id << "\n";

    file.close();
}

std::string Data::format_double(double value, int precision) const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(precision) << value;

    return oss.str();
}