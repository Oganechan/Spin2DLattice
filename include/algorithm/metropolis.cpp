#include "metropolis.hpp"
#include <cstdint>

algorithm::Metropolis::Metropolis(lattice::Atoms &atoms,
                                  const physics::Calculator &calculator)
    : atoms_(atoms), calculator_(calculator) {
    Random::uniform_real<double>();
}

void algorithm::Metropolis::step() {
    const int32_t atom_id = atoms_.select_random_magnetic_atom();
    const std::array<double, 3> new_spin = atoms_.small_rotate_spin(atom_id);

    const double energy_diff =
        calculator_.calculate_flip_energy_difference(atom_id, new_spin);

    if (accept_change(energy_diff))
        atoms_.set_spin(atom_id, new_spin);
}

void algorithm::Metropolis::sweep() {
    int32_t step_count = atoms_.get_magnetic_count();

    while (step_count-- > 0)
        step();
}

void algorithm::Metropolis::sweep(int32_t step_count) {
    while (step_count-- > 0)
        step();
}

bool algorithm::Metropolis::accept_change(double energy_difference) const {
    if (energy_difference <= 0)
        return true;

    double random_value = Random::uniform_real<double>();
    double threshold =
        std::exp(-energy_difference / calculator_.get_temperature());

    return random_value <= threshold;
}