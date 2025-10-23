#include "calculator.h"

physics::Calculator::Calculator(const lattice::Atoms &atoms)
    : atoms_(atoms),
      geometry_(atoms.get_geometry()),
      neighbor_table_(geometry_.get_neighbor_table()),
      atom_count_(geometry_.get_atom_count()),
      shell_count_(geometry_.get_shell_count()),
      exchange_constants_(atoms_.get_config().get<std::vector<double>>("physical.exchange_constants")),
      external_magnetic_field_(atoms.get_config().get<std::array<double, 3>>("physical.external_magnetic_field")) {}

// === ENERGY CALCULATIONS ===

double physics::Calculator::calculate_total_energy() const
{
    double total_energy = 0.0;

    auto magnetic_atoms = atoms_.get_magnetic_atoms();
    for (int32_t atom_id : magnetic_atoms)
        total_energy += calculate_atom_energy(atom_id);

    return total_energy / 2.0;
}

double physics::Calculator::calculate_atom_energy(int32_t atom_id) const
{
    if (!atoms_.get_magnetic_state(atom_id))
        return 0.0;

    const auto &atom_spin = atoms_.get_spin(atom_id);
    double energy = -calculate_spin_dot_product(atom_spin, external_magnetic_field_);

    for (int32_t shell_index = 0; shell_index < shell_count_; ++shell_index)
    {
        double J = get_exchange_constant(shell_index);
        for (int32_t neighbor_id : neighbor_table_[atom_id][shell_index])
            if (atoms_.get_magnetic_state(neighbor_id))
                energy += J * calculate_spin_dot_product(atom_spin, atoms_.get_spin(neighbor_id));
    }

    return energy;
}

double physics::Calculator::calculate_pair_energy(int32_t first_atom_id, int32_t second_atom_id) const
{
    if (!atoms_.get_magnetic_state(first_atom_id) || !atoms_.get_magnetic_state(second_atom_id))
        return 0.0;

    const auto &spin1 = atoms_.get_spin(first_atom_id);
    const auto &spin2 = atoms_.get_spin(second_atom_id);

    double distance = geometry_.get_distance(first_atom_id, second_atom_id);
    double J = get_exchange_constant(static_cast<int32_t>(distance));

    return J * calculate_spin_dot_product(spin1, spin2);
}

double physics::Calculator::calculate_flip_energy_difference(int32_t atom_id) const
{
    if (!atoms_.get_magnetic_state(atom_id))
        return 0.0;

    double original_energy = calculate_atom_energy(atom_id);

    double x = Random::uniform_real(-1.0, 1.0);
    double y = Random::uniform_real(-1.0, 1.0);
    double z = Random::uniform_real(-1.0, 1.0);
    double norm = std::sqrt(x * x + y * y + z * z);
    std::array<double, 3> new_spin = {x / norm, y / norm, z / norm};

    double new_energy = -calculate_spin_dot_product(new_spin, external_magnetic_field_);

    for (int32_t shell_index = 0; shell_index < shell_count_; ++shell_index)
    {
        double J = get_exchange_constant(shell_index);
        for (int32_t neighbor_id : neighbor_table_[atom_id][shell_index])
            if (atoms_.get_magnetic_state(neighbor_id))
                new_energy += J * calculate_spin_dot_product(new_spin, atoms_.get_spin(neighbor_id));
    }

    return new_energy - original_energy;
}

// === MAGNETIZATION CALCULATIONS ===

double physics::Calculator::calculate_total_magnetization() const
{
    int32_t magnetic_count = atoms_.get_magnetic_count();
    if (magnetic_count == 0)
        return 0.0;

    std::array<double, 3> total_magnetization = {0.0, 0.0, 0.0};

    auto magnetic_atoms = atoms_.get_magnetic_atoms();
    for (int32_t atom_id : magnetic_atoms)
    {
        const auto &spin = atoms_.get_spin(atom_id);
        total_magnetization[0] += spin[0];
        total_magnetization[1] += spin[1];
        total_magnetization[2] += spin[2];
    }

    double magnitude = std::sqrt(total_magnetization[0] * total_magnetization[0] +
                                 total_magnetization[1] * total_magnetization[1] +
                                 total_magnetization[2] * total_magnetization[2]);

    return magnitude / static_cast<double>(magnetic_count);
}