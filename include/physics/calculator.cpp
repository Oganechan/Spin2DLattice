#include "calculator.hpp"
#include "../core/config.hpp"
#include <cstdint>
#include <vector>

physics::Calculator::Calculator(const lattice::Atoms &atoms)
    : atoms_(atoms), geometry_(atoms.get_geometry()),
      neighbor_table_(geometry_.get_neighbor_table()),
      atom_count_(geometry_.get_atom_count()),
      shell_count_(geometry_.get_shell_count()),
      exchange_constants_(atoms_.get_config().get<std::vector<double>>(
          "physical.exchange_constants")),
      external_magnetic_field_(atoms.get_config().get<std::array<double, 3>>(
          "physical.external_magnetic_field")) {
    magnetic_atoms_.reserve(atom_count_);
}

// === ENERGY CALCULATIONS ===

double physics::Calculator::calculate_total_energy() const {
    double total_energy = 0.0;
    auto magnetic_atoms = atoms_.get_magnetic_atoms();

    for (int32_t atom_id : magnetic_atoms) {
        double energy = 0;
        const auto &atom_spin = atoms_.get_spin(atom_id);

        energy -=
            calculate_spin_dot_product(atom_spin, external_magnetic_field_);

        for (int shell_index = 0; shell_index < shell_count_; ++shell_index) {
            const double J = exchange_constants_[shell_index];
            const auto &neighbors = neighbor_table_[atom_id][shell_index];

            for (int32_t neighbor_id : neighbors)
                if (atoms_.get_magnetic_state(neighbor_id))
                    energy -= J * calculate_spin_dot_product(
                                      atom_spin, atoms_.get_spin(neighbor_id));
        }
        total_energy += energy;
    }

    return total_energy / 2.0;
}

double physics::Calculator::calculate_atom_energy(int32_t atom_id) const {
    if (!atoms_.get_magnetic_state(atom_id))
        return 0.0;

    const auto &atom_spin = atoms_.get_spin(atom_id);
    double energy =
        -calculate_spin_dot_product(atom_spin, external_magnetic_field_);

    for (int32_t shell_index = 0; shell_index < shell_count_; ++shell_index) {
        const double J = get_exchange_constant(shell_index);
        const auto &neighbors = neighbor_table_[atom_id][shell_index];

        for (int32_t neighbor_id : neighbors)
            if (atoms_.get_magnetic_state(neighbor_id))
                energy -= J * calculate_spin_dot_product(
                                  atom_spin, atoms_.get_spin(neighbor_id));
    }

    return energy;
}

double
physics::Calculator::calculate_pair_energy(int32_t first_atom_id,
                                           int32_t second_atom_id) const {
    if (!atoms_.get_magnetic_state(first_atom_id) ||
        !atoms_.get_magnetic_state(second_atom_id))
        return 0.0;

    const auto &spin1 = atoms_.get_spin(first_atom_id);
    const auto &spin2 = atoms_.get_spin(second_atom_id);

    double distance = geometry_.get_distance(first_atom_id, second_atom_id);
    double J = get_exchange_constant(static_cast<int32_t>(distance));

    return -J * calculate_spin_dot_product(spin1, spin2);
}

double physics::Calculator::calculate_flip_energy_difference(
    int32_t atom_id, std::array<double, 3> new_spin) const {
    if (!atoms_.get_magnetic_state(atom_id))
        return 0.0;

    const auto &old_spin = atoms_.get_spin(atom_id);

    const double spin_diff_x = new_spin[0] - old_spin[0];
    const double spin_diff_y = new_spin[1] - old_spin[1];
    const double spin_diff_z = new_spin[2] - old_spin[2];

    double energy_diff = -(spin_diff_x * external_magnetic_field_[0] +
                           spin_diff_y * external_magnetic_field_[1] +
                           spin_diff_z * external_magnetic_field_[2]);

    for (int32_t shell_index = 0; shell_index < shell_count_; ++shell_index) {
        const double J = exchange_constants_[shell_index];
        const auto &neighbors = neighbor_table_[atom_id][shell_index];

        for (int32_t neighbor_id : neighbors)
            if (atoms_.get_magnetic_state(neighbor_id)) {
                const auto &neighbor_spin = atoms_.get_spin(neighbor_id);
                energy_diff -= J * (spin_diff_x * neighbor_spin[0] +
                                    spin_diff_y * neighbor_spin[1] +
                                    spin_diff_z * neighbor_spin[2]);
            }
    }

    return energy_diff;
}

// === MAGNETIZATION CALCULATIONS ===

double physics::Calculator::calculate_total_magnetization() const {
    const auto &magnetic_atoms = get_magnetic_atoms_cached();
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

const std::vector<int32_t> &
physics::Calculator::get_magnetic_atoms_cached() const {
    if (!magnetic_atoms_valid_) {
        update_magnetic_cache();
    }
    return magnetic_atoms_;
}

void physics::Calculator::update_magnetic_cache() const {
    magnetic_atoms_.clear();
    const auto &magnetic_mask = atoms_.get_magnetic_mask();

    for (int32_t i = 0; i < atom_count_; ++i) {
        if (magnetic_mask[i]) {
            magnetic_atoms_.push_back(i);
        }
    }
    magnetic_atoms_valid_ = true;
}