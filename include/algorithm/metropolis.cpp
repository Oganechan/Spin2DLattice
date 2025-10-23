#include "metropolis.h"

algorithm::Metropolis::Metropolis(lattice::Atoms &atoms, double temperature)
    : atoms_(atoms), calculator_(atoms), temperature_(temperature), beta_(1.0 / temperature) {}

void algorithm::Metropolis::step()
{
    int32_t atom_id = select_random_magnetic_atom();
    std::array<double, 3> new_spin = atoms_.generate_random_spin();
    double energy_diff = calculator_.calculate_flip_energy_difference(atom_id, new_spin);

    if (accept_change(energy_diff))
        atoms_.set_spin(atom_id, new_spin);
}

void algorithm::Metropolis::sweep()
{
    const int32_t magnetic_count = atoms_.get_magnetic_count();

    for (int32_t i = 0; i < magnetic_count; ++i)
        step();
}

int32_t algorithm::Metropolis::select_random_magnetic_atom() const
{
    const auto &magnetic_atoms = atoms_.get_magnetic_atoms();
    return magnetic_atoms[Random::uniform_int<int32_t>(0, magnetic_atoms.size() - 1)];
}