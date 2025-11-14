#include "metropolis.h"

algorithm::Metropolis::Metropolis(lattice::Atoms &atoms)
    : atoms_(atoms), calculator_(atoms),
      temperature_(atoms_.get_config().get<double>("physical.temperature")),
      beta_(1.0 / temperature_) {}

void algorithm::Metropolis::step() {
    int32_t atom_id = atoms_.select_random_magnetic_atom();
    std::array<double, 3> new_spin = atoms_.generate_random_spin();
    double energy_diff =
        calculator_.calculate_flip_energy_difference(atom_id, new_spin);

    if (accept_change(energy_diff))
        atoms_.set_spin(atom_id, new_spin);
}

void algorithm::Metropolis::sweep() {
    int32_t magnetic_count = atoms_.get_magnetic_count();

    while (magnetic_count-- > 0)
        step();
}