#pragma once

#include <variant>
#include "../geometry.h"
#include "../../utils/random.h"

namespace lattice
{

    using SpinVariant = std::variant<int32_t, std::array<double, 3>>;

    class SpinModelBase
    {
    public:
        virtual ~SpinModelBase() = default;

        // === MODEL METADATA ===
        virtual SpinModel get_type() const = 0;

        // === DEFECT MANAGEMENT ===
        virtual void set_magnetic(int32_t atom_index, bool magnetic) = 0;
        virtual bool is_magnetic(int32_t atom_index) const = 0;
        virtual const std::vector<bool> &get_magnetic_mask() const = 0;
        virtual void set_random_defects(double concentration) = 0;

        // === SYSTEM STATISTICS ===
        virtual int32_t count_magnetic_atoms() const = 0;
        virtual int32_t count_defects() const = 0;
        virtual std::vector<int32_t> get_magnetic_atom_indices() const = 0;
        virtual std::vector<int32_t> get_defect_indices() const = 0;

        // === WORKING WITH SPINS ===
        virtual SpinVariant get_spin(int32_t atom_index) const = 0;
        virtual void set_spin(int32_t atom_index, const SpinVariant &value) = 0;
        virtual void flip_spins(const std::vector<int32_t> &indices) = 0;

        // === INITIALIZING CONFIGURATIONS ===
        virtual void random_initialize() = 0;
    };

} // namespace lattice