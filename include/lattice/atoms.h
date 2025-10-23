#pragma once

#include "geometry.h"
#include "../utils/random.h"

namespace lattice
{

    class Atoms
    {
    public:
        explicit Atoms(const Config &config);

        // === DEFECT MANAGEMENT ===

        void set_magnetic(int32_t atom_index, bool magnetic);
        bool is_magnetic(int32_t atom_index) const;
        const std::vector<bool> get_magnetic_mask();
        void set_random_defects(double concentration);

        // === SYSTEM STATISTICS ===

        int32_t count_magnetic_atoms() const;
        int32_t count_defects() const;
        std::vector<int32_t> get_magnetic_atom_indices() const;
        std::vector<int32_t> get_defect_indices() const;

        // === WORKING WITH SPINS ===

        std::array<double, 3> get_spin(int32_t atom_index) const;
        void set_spin(int32_t atom_index, const std::array<double, 3> &value);
        void flip_spins(const std::vector<int32_t> &indices);

        // === INITIALIZING CONFIGURATIONS ===

        void random_initialize();
        void ferromagnetic_initialize();
        void antiferromagnetic_initialize();

        const Config &config() const { return config_; }
        const Geometry &geometry() const { return *geometry_; }
        int32_t get_num_atoms() const { return geometry_->get_num_atoms(); }

    private:
        const Config &config_;
        std::unique_ptr<Geometry> geometry_;

        std::vector<std::array<double, 3>> spins_;
        std::vector<bool> magnetic_mask_;

        void initialize_geometry();
        void initialize_spins();
    };

} // namespace lattice