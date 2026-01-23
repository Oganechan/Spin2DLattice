#ifndef DATA_HPP
#define DATA_HPP

#include "../physics/calculator.hpp"
#include <string>
#include <vector>

class Data {
  public:
    explicit Data(const physics::Calculator &calculator,
                  const std::string &output_path);

    void measure();
    void save();
    void reset();

    inline double get_mean_energy() const { return mean_energy_; }
    inline double get_mean_energy_sq() const { return mean_energy_sq; }

    inline double get_mean_order_parameter() const {
        return mean_order_parameter_;
    }
    inline double get_mean_order_parameter_sq() const {
        return mean_order_parameter_sq_;
    }
    inline double get_mean_order_parameter_4th() const {
        return mean_order_parameter_4th_;
    }

    inline double get_specific_heat() const { return specific_heat_; }

    inline double get_order_susceptibility() const {
        return order_susceptibility_;
    }

    inline double get_order_binder() const { return order_binder_; }

  private:
    const physics::Calculator &calculator_;
    const std::string &output_path_;

    std::vector<double> energies_;
    std::vector<double> magnetizations_;
    std::vector<double> order_parameters_;

    double mean_energy_ = 0.0;
    double mean_energy_sq = 0.0;

    double mean_order_parameter_ = 0.0;     // ⟨φ⟩
    double mean_order_parameter_abs_ = 0.0; // ⟨|φ|⟩
    double mean_order_parameter_sq_ = 0.0;  // ⟨φ²⟩
    double mean_order_parameter_4th_ = 0.0; // ⟨φ⁴⟩

    double specific_heat_ = 0.0;
    double order_susceptibility_ = 0.0;
    double order_binder_ = 0.0;

    void compute_statistics();
};

#endif // DATA_HPP