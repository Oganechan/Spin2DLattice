#pragma once

#include "../utils/random.h"
#include "geometry.h"

class Atoms
{
public:
    explicit Atoms(Lattice &lattice)
        : num_atoms_(lattice.get_num_atoms())
    {
        spins_.reserve(num_atoms_);
        magnetic_indices_.reserve(num_atoms_);
    }

    void set_spin(int32_t index, int32_t spin);
    void set_paramagnetic(int32_t index);
    void set_magnetic(int32_t index);

    const std::vector<int32_t> &get_spins() const { return spins_; }
    const std::vector<int32_t> &get_magnetic_atoms() const { return magnetic_indices_; }
    int32_t get_spin(int32_t index) const { return spins_[index]; }
    int32_t get_num_magnetic() const { return num_magnetic_atoms; }

    void random_initialize(double concentration);

    void flip_spin(int32_t index);
    bool is_magnetic(int32_t index) const;

private:
    const int32_t num_atoms_;
    int32_t num_magnetic_atoms;

    std::vector<int32_t> spins_, magnetic_indices_;
};