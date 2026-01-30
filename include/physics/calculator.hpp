#ifndef CALCULATOR_HPP
#define CALCULATOR_HPP

#include "../lattice/atoms.hpp"
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace physics {

class Calculator {
  public:
    explicit Calculator(const lattice::Atoms &atoms, const Config &config);

    double calculate_total_energy() const;
    double calculate_atom_energy(int32_t atom_id) const;
    double calculate_flip_energy_difference(
        int32_t atom_id, const std::array<double, 3> &new_spin) const;
    double calculate_flip_energy_difference(
        int32_t atom_id,
        const std::unique_ptr<lattice::BaseSpin> &new_spin) const;
    double calculate_total_magnetization() const;
    double calculate_total_order() const;

    inline void set_temperature(double temperature) {
        temperature_ = temperature;
    }

    inline void set_concentration(double concentration) {
        concentration_ = concentration;
    }

    inline double get_temperature() const { return temperature_; }
    inline double get_concentration() const { return concentration_; }
    inline double get_exchange_constant(int32_t shell_index) const {
        return exchange_constants_[shell_index];
    }
    inline const lattice::Atoms &get_atoms() const { return atoms_; }

  private:
    const lattice::Atoms &atoms_;
    const lattice::Geometry &geometry_;

    const std::array<double, 3> external_magnetic_field_;
    const std::vector<double> exchange_constants_;
    const double anisotropy_constant_;

    double temperature_;
    double concentration_;

    double dot_product(const lattice::BaseSpin &spin1,
                       const lattice::BaseSpin &spin2) const;
    double dot_product(const lattice::BaseSpin &spin1,
                       const std::array<double, 3> &spin2) const;
    double calculate_anisotropy_energy(const lattice::BaseSpin &spin) const;
};

} // namespace physics

#endif // CALCULATOR_HPP