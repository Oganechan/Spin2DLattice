#include "../../include/lattice/atoms.h"

void Atoms::set_spin(int32_t index, int32_t spin)
{
    if (index < 0 || index >= num_atoms_)
        throw std::invalid_argument("Index is out of range num_atoms");

    if (spin < -1 || spin > 1)
        throw std::invalid_argument("Spin can be only -1, 0 or 1");

    spins_[index] = spin;
}

void Atoms::set_paramagnetic(int32_t index)
{
    if (index < 0 || index >= num_atoms_)
        throw std::invalid_argument("Index is out of range num_atoms");

    spins_[index] = 0;
}

void Atoms::set_magnetic(int32_t index)
{
    if (index < 0 || index >= num_atoms_)
        throw std::invalid_argument("Index is out of range num_atoms");

    std::bernoulli_distribution dist(0.5);
    spins_[index] = dist(Random::get_rng()) ? 1 : -1;
}

void Atoms::random_initialize(double concentration)
{
    if (concentration < 0 || concentration > 1)
        throw std::invalid_argument("Incorrect concentration value");

    spins_.clear();
    magnetic_indices_.clear();

    std::bernoulli_distribution dist(0.5);
    for (int32_t i = 0; i < num_atoms_; ++i)
    {
        spins_.push_back(dist(Random::get_rng()) ? 1 : -1);
        magnetic_indices_.push_back(i);
    }

    int32_t zero_spins_count = static_cast<int32_t>((1.0 - concentration) * static_cast<double>(num_atoms_));
    if (zero_spins_count == 0 || magnetic_indices_.empty())
        return;
    for (int32_t i = 0; i < zero_spins_count && !magnetic_indices_.empty(); ++i)
    {
        std::uniform_int_distribution<int32_t> dist(0, magnetic_indices_.size() - 1);
        int32_t random_index = dist(Random::get_rng());
        spins_[magnetic_indices_[random_index]] = 0;
        magnetic_indices_.erase(magnetic_indices_.begin() + random_index);
    }
}

void Atoms::flip_spin(int32_t index)
{
    if (index < 0 || index >= num_atoms_)
        throw std::invalid_argument("Index is out of range num_atoms");

    if (!is_magnetic(index))
        return;

    spins_[index] *= -1;
}

bool Atoms::is_magnetic(int32_t index) const
{
    if (index < 0 || index >= num_atoms_)
        throw std::invalid_argument("Index is out of range num_atoms");

    if (spins_[index] == 0)
        return false;
    else
        return true;
}