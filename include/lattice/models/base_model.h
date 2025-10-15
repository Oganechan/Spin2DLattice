#pragma once

#include "../geometry.h"
#include "../../utils/random.h"

namespace lattice
{

    class SpinModelBase
    {
    public:
        virtual ~SpinModelBase() = default;

        // Set magnetic state for specific atom (true = magnetic, false = paramagnetic defects)
        virtual void set_magnetic(int32_t atom_index, bool magnetic) = 0;

        // Check if atom is magnetic
        virtual bool is_magnetic(int32_t atom_index) const = 0;

        // Get reference to the magnetic mask vector for bulk operations
        virtual const std::vector<bool> &get_magnetic_mask() const = 0;

        // Create random paramagnetic defects with given concentration [0,1]
        virtual void set_random_defects(double concentration) = 0;

        // Randomly initialize spins
        virtual void random_initialize() = 0;

        // Initialize all magnetic spins in ferromagnetic alignment
        virtual void ferromagnetic_initialize() = 0;

        // Initialize all magnetic spins in antiferromagnetic pattern
        virtual void antiferromagnetic_initialize() = 0;

        // Get the type of spin model (Ising, Heisenberg, etc.)
        virtual SpinModel get_type() const = 0;
    };

}
