#include "calculator.hpp"
#include <cstdint>
#include <vector>

physics::Calculator::Calculator(const lattice::Atoms &atoms,
                                const Config &config)
    : atoms_(atoms), geometry_(atoms.get_geometry()),
      external_magnetic_field_(config.get<std::array<double, 3>>(
          "physical.external_magnetic_field")),
      exchange_constants_(
          config.get<std::vector<double>>("physical.exchange_constants")),
      temperature_(config.get<double>("physical.temperature")),
      concentration_(config.get<double>("physical.concentration")) {}

// Calculates the total energy of the spin system
double physics::Calculator::calculate_total_energy() const {
    double total_energy = 0.0;

    for (int32_t atom_id : atoms_.get_magnetic_atoms()) {
        const auto &spin = atoms_.get_spin(atom_id);
        total_energy -=
            calculate_spin_dot_product(spin, external_magnetic_field_);
    }

    for (int32_t atom_id : atoms_.get_magnetic_atoms()) {
        const auto &spin1 = atoms_.get_spin(atom_id);

        for (int32_t shell_index = 0; shell_index < geometry_.get_shell_count();
             ++shell_index) {
            const double J = get_exchange_constant(shell_index);

            for (int32_t neighbor_id :
                 geometry_.get_neighbor_table()[atom_id][shell_index]) {
                if (atoms_.get_magnetic_state(neighbor_id) &&
                    neighbor_id > atom_id) {
                    const auto &spin2 = atoms_.get_spin(neighbor_id);
                    total_energy -=
                        J * calculate_spin_dot_product(spin1, spin2);
                }
            }
        }
    }

    return total_energy;
}

// Calculates the local energy contribution for a specific atom
double physics::Calculator::calculate_atom_energy(int32_t atom_id) const {
    if (!atoms_.get_magnetic_state(atom_id))
        return 0.0;

    const auto &atom_spin = atoms_.get_spin(atom_id);
    double energy =
        -calculate_spin_dot_product(atom_spin, external_magnetic_field_);

    for (int32_t shell_index = 0; shell_index < geometry_.get_shell_count();
         ++shell_index) {
        const double J = get_exchange_constant(shell_index);
        const auto &neighbors =
            geometry_.get_neighbor_table()[atom_id][shell_index];

        for (int32_t neighbor_id : neighbors)
            if (atoms_.get_magnetic_state(neighbor_id))
                energy -= J * calculate_spin_dot_product(
                                  atom_spin, atoms_.get_spin(neighbor_id));
    }

    return energy;
}

// Calculates the energy difference if the specified atom's spin were
// flipped
double physics::Calculator::calculate_flip_energy_difference(
    int32_t atom_id, const std::array<double, 3> &new_spin) const {
    if (!atoms_.get_magnetic_state(atom_id))
        return 0.0;

    const auto &old_spin = atoms_.get_spin(atom_id);

    const double spin_diff_x = new_spin[0] - old_spin[0];
    const double spin_diff_y = new_spin[1] - old_spin[1];
    const double spin_diff_z = new_spin[2] - old_spin[2];

    double energy_diff = -(spin_diff_x * external_magnetic_field_[0] +
                           spin_diff_y * external_magnetic_field_[1] +
                           spin_diff_z * external_magnetic_field_[2]);

    for (int32_t shell_index = 0; shell_index < geometry_.get_shell_count();
         ++shell_index) {
        const double J = exchange_constants_[shell_index];

        for (int32_t neighbor_id :
             geometry_.get_neighbor_table()[atom_id][shell_index])
            if (atoms_.get_magnetic_state(neighbor_id)) {
                const auto &neighbor_spin = atoms_.get_spin(neighbor_id);
                energy_diff -= J * (spin_diff_x * neighbor_spin[0] +
                                    spin_diff_y * neighbor_spin[1] +
                                    spin_diff_z * neighbor_spin[2]);
            }
    }

    return energy_diff;
}

// Calculates the total magnetization vector magnitude of the system
double physics::Calculator::calculate_total_magnetization() const {
    const auto &magnetic_atoms = atoms_.get_magnetic_atoms();
    int32_t magnetic_count = atoms_.get_magnetic_count();

    if (magnetic_count == 0)
        return 0.0;

    double total_x = 0.0, total_y = 0.0, total_z = 0.0;

    for (int32_t atom_id : magnetic_atoms) {
        const auto &spin = atoms_.get_spin(atom_id);
        total_x += spin[0];
        total_y += spin[1];
        total_z += spin[2];
    }

    return std::sqrt(total_x * total_x + total_y * total_y +
                     total_z * total_z) /
           magnetic_count;
}
