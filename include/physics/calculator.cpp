#include "calculator.hpp"
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

physics::Calculator::Calculator(const lattice::Atoms &atoms,
                                const Config &config)
    : atoms_(atoms), geometry_(atoms.get_geometry()),
      external_magnetic_field_(config.get<std::array<double, 3>>(
          "physical.external_magnetic_field")),
      exchange_constants_(
          config.get<std::vector<double>>("physical.exchange_constants")),
      anisotropy_constant_(config.get<double>("physical.anisotropy_constant")),
      temperature_(config.get<double>("physical.temperature")),
      concentration_(config.get<double>("physical.concentration")) {}

// Calculates the total energy of the spin system
double physics::Calculator::calculate_total_energy() const {
    double total_energy = 0.0;

    for (int32_t atom_id : atoms_.get_magnetic_atoms()) {
        const auto &spin = atoms_.get_spin(atom_id);
        total_energy +=
            -calculate_spin_dot_product(spin, external_magnetic_field_) -
            anisotropy_constant_ * spin[2] * spin[2];
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
                    total_energy +=
                        -J * calculate_spin_dot_product(spin1, spin2);
                }
            }
        }
    }

    return total_energy / atoms_.get_geometry().get_atom_count();
}

// Calculates the local energy contribution for a specific atom
double physics::Calculator::calculate_atom_energy(int32_t atom_id) const {
    if (!atoms_.get_magnetic_state(atom_id))
        return 0.0;

    const auto &atom_spin = atoms_.get_spin(atom_id);
    double energy =
        -calculate_spin_dot_product(atom_spin, external_magnetic_field_) -
        anisotropy_constant_ * atom_spin[2] * atom_spin[2];

    for (int32_t shell_index = 0; shell_index < geometry_.get_shell_count();
         ++shell_index) {
        const double J = get_exchange_constant(shell_index);
        const auto &neighbors =
            geometry_.get_neighbor_table()[atom_id][shell_index];

        for (int32_t neighbor_id : neighbors)
            if (atoms_.get_magnetic_state(neighbor_id))
                energy += -J * calculate_spin_dot_product(
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

    double energy_diff = -calculate_spin_dot_product(
        {spin_diff_x, spin_diff_y, spin_diff_z}, external_magnetic_field_);

    energy_diff -= anisotropy_constant_ *
                   (new_spin[2] * new_spin[2] - old_spin[2] * old_spin[2]);

    for (int32_t shell_index = 0; shell_index < geometry_.get_shell_count();
         ++shell_index) {
        const double J = get_exchange_constant(shell_index);

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
    if (atoms_.get_magnetic_count() == 0)
        return 0.0;

    const auto &magnetic_atoms = atoms_.get_magnetic_atoms();
    double total_x = 0.0, total_y = 0.0, total_z = 0.0;

    for (int32_t atom_id : magnetic_atoms) {
        const auto &spin = atoms_.get_spin(atom_id);
        total_x += spin[0];
        total_y += spin[1];
        total_z += spin[2];
    }

    return std::sqrt(total_x * total_x + total_y * total_y +
                     total_z * total_z) /
           geometry_.get_atom_count();
}

double physics::Calculator::calculate_total_order() const {
    if (atoms_.get_magnetic_count() == 0)
        return 0.0;

    int32_t K = geometry_.get_sublattice_count();
    int32_t N_total = geometry_.get_atom_count();

    std::vector<std::array<double, 3>> sublattice_sum(K, {0.0, 0.0, 0.0});
    std::vector<int32_t> sublattice_count(K, 0);

    for (int32_t atom_id : atoms_.get_magnetic_atoms()) {
        const auto &spin = atoms_.get_spin(atom_id);
        int32_t sublattice = geometry_.get_atom_sublattices()[atom_id];

        sublattice_sum[sublattice][0] += spin[0];
        sublattice_sum[sublattice][1] += spin[1];
        sublattice_sum[sublattice][2] += spin[2];
        sublattice_count[sublattice]++;
    }

    std::vector<std::array<double, 3>> sublattice_magnetization(K);
    for (int r = 0; r < K; ++r) {
        if (sublattice_count[r] > 0) {
            sublattice_magnetization[r][0] =
                sublattice_sum[r][0] / sublattice_count[r];
            sublattice_magnetization[r][1] =
                sublattice_sum[r][1] / sublattice_count[r];
            sublattice_magnetization[r][2] =
                sublattice_sum[r][2] / sublattice_count[r];
        }
    }

    double sum_sq = 0.0;
    for (int r = 0; r < K; ++r) {
        double mx = sublattice_magnetization[r][0];
        double my = sublattice_magnetization[r][1];
        double mz = sublattice_magnetization[r][2];
        sum_sq += mx * mx + my * my + mz * mz;
    }

    return std::sqrt(sum_sq / K);
}
