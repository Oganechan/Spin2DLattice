#pragma once

#include "../lattice/atoms.h"
#include "../core/config.h"

namespace physics
{

    class Calculator
    {
    public:
        explicit Calculator(const lattice::Atoms &atoms);

        // === ENERGY CALCULATIONS ===

        // Calculates the total energy of the spin system
        double calculate_total_energy() const;

        // Calculates the local energy contribution for a specific atom
        double calculate_atom_energy(int32_t atom_id) const;

        // Calculates the exchange energy between two specific atoms
        double calculate_pair_energy(int32_t first_atom_id, int32_t second_atom_id) const;

        // Calculates the energy difference if the specified atom's spin were flipped
        double calculate_flip_energy_difference(int32_t atom_id, std::array<double, 3> new_spin) const;

        // === MAGNETIZATION CALCULATIONS ===

        // Calculates the total magnetization vector magnitude of the system
        double calculate_total_magnetization() const;

    private:
        const lattice::Atoms &atoms_;
        const lattice::Geometry &geometry_;
        const std::vector<std::vector<std::vector<int32_t>>> &neighbor_table_;

        const int32_t atom_count_;
        const int32_t shell_count_;
        const std::vector<double> exchange_constants_;
        const std::array<double, 3> external_magnetic_field_;

        inline double get_exchange_constant(int32_t shell_index) const
        {
            if (shell_index < 0 || shell_index >= exchange_constants_.size())
                throw std::out_of_range("Shell index out of range");

            return exchange_constants_[shell_index];
        }

        inline double calculate_spin_dot_product(const std::array<double, 3> &spin1,
                                                 const std::array<double, 3> &spin2) const
        {
            return spin1[0] * spin2[0] + spin1[1] * spin2[1] + spin1[2] * spin2[2];
        }
    };

} // namespace physics