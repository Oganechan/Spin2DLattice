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

        // Sets the magnetic state of an atom (true = magnetic, false = non-magnetic)
        void set_magnetic_state(int32_t atom_id, bool is_magnetic);

        // Returns the magnetic state of the specified atom
        bool get_magnetic_state(int32_t atom_id) const;

        // Returns a vector representing the magnetic mask for all atoms
        const std::vector<bool> get_magnetic_mask();

        // Randomly sets atoms as non-magnetic atoms with given concentration
        void set_random_defects(double defect_concentration);

        // === SYSTEM STATISTICS ===

        // Returns the total number of magnetic atoms in the system
        int32_t get_magnetic_count() const;

        // Returns the total number of non-magnetic atoms in the system
        int32_t get_defect_count() const;

        // Returns a list of atom IDs that are magnetic
        std::vector<int32_t> get_magnetic_atoms() const;

        // Returns a list of atom IDs that are non-magnetic
        std::vector<int32_t> get_defect_atoms() const;

        // === WORKING WITH SPINS ===

        // Returns the spin vector of the specified atom
        std::array<double, 3> get_spin(int32_t atom_id) const;

        // Sets the spin vector for the specified atom
        void set_spin(int32_t atom_id, const std::array<double, 3> &spin_vector);

        // Flips the spin direction for the specified atoms
        void flip_spins(const std::vector<int32_t> &atom_ids);

        // === INITIALIZING CONFIGURATIONS ===

        // Initializes all magnetic spins with random orientations
        void initialize_random();

        // Initializes all magnetic spins in ferromagnetic alignment (+z direction)
        void initialize_ferromagnetic();

        // Initializes magnetic spins in antiferromagnetic pattern (+-z alternating direction)
        void initialize_antiferromagnetic();

        // Returns the configuration object
        const Config &get_config() const { return config_; }

        // Returns the geometry object
        const Geometry &get_geometry() const { return geometry_; }

    private:
        const Config &config_;
        const Geometry geometry_;

        std::vector<std::array<double, 3>> spin_vectors_;
        std::vector<bool> magnetic_mask_;

        void initialize_geometry();
        void initialize_spins();
    };

} // namespace lattice