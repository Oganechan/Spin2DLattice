#ifndef DATA_HPP
#define DATA_HPP

#include "../physics/calculator.hpp"
#include <string>

class Data {
  public:
    Data(const physics::Calculator &calculator, const std::string &output_path);

    void measure();
    void save();
    void reset();

    inline double get_mean_energy() const { return mean_energy_; }
    inline double get_mean_magnetization() const { return mean_magnetization_; }

  private:
    const physics::Calculator &calculator_;

    const std::string &output_path_;

    std::vector<double> energies_;
    std::vector<double> magnetizations_;

    double mean_energy_ = 0.0;
    double mean_magnetization_ = 0.0;

    void compute_statistics();
};

#endif // DATA_HPP