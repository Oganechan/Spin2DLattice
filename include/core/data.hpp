#ifndef DATA_HPP
#define DATA_HPP

#include "../physics/calculator.hpp"
#include "config.hpp"
#include "measurements.hpp"
#include <cstdint>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

class Data {
  public:
    explicit Data(const physics::Calculator &calculator, const Config &config,
                  const std::string &base_output_dir);

    void measure();
    void save_statistics();
    void save_finale();

  private:
    const physics::Calculator &calculator_;
    const Config &config_;

    fs::path output_dir_;
    std::string simulation_name_;

    Measurements measurements_;

    std::string generate_simulation_name() const;

    void initialize();

    void create_directory_structure();
    void save_config() const;
    void save_geometry() const;
    void save_neighbor_table() const;

    void reset();
    void compute_statistics();

    std::string format_double(double value, int32_t precision = 8) const;
};

#endif // DATA_HPP