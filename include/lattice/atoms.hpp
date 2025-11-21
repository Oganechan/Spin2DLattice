#ifndef ATOMS_HPP
#define ATOMS_HPP

#include "../utils/random.hpp"
#include "geometry.hpp"
#include <cstdint>
#include <vector>

namespace lattice {

class Atoms {
  public:
    explicit Atoms(const Config &config);

    // === DEFECT MANAGEMENT ===

    // Sets the magnetic state of an atom (true = magnetic, false =
    // non-magnetic)
    inline void set_magnetic_state(int32_t atom_id, bool is_magnetic) {
        magnetic_mask_[atom_id] = is_magnetic;
    }

    // Returns the magnetic state of the specified atom
    inline bool get_magnetic_state(int32_t atom_id) const {
        return magnetic_mask_[atom_id];
    }

    // Returns a vector representing the magnetic mask for all atoms
    const std::vector<bool> &get_magnetic_mask() const {
        return magnetic_mask_;
    }

    // Randomly sets atoms as non-magnetic atoms with given concentration
    void set_random_defects(double defect_concentration);

    // === SYSTEM STATISTICS ===

    // Returns the total number of magnetic atoms in the system
    inline int32_t get_magnetic_count() const { return magnetic_count_; }

    // Returns the total number of non-magnetic atoms in the system
    inline int32_t get_defect_count() const {
        return geometry_.get_atom_count() - magnetic_count_;
    }

    // Returns a list of atom IDs that are magnetic
    inline const std::vector<int32_t> &get_magnetic_atoms() const {
        return magnetic_atoms_;
    }

    // Returns a list of atom IDs that are non-magnetic
    inline const std::vector<int32_t> &get_defect_atoms() const {
        return defect_atoms_;
    }

    // === WORKING WITH SPINS ===

    // Returns the spin vector of the specified atom
    inline std::array<double, 3> get_spin(int32_t atom_id) const {
        return spin_vectors_[atom_id];
    }

    // Sets the spin vector for the specified atom
    inline void set_spin(int32_t atom_id,
                         const std::array<double, 3> &spin_vector) {
        double norm = std::sqrt(spin_vector[0] * spin_vector[0] +
                                spin_vector[1] * spin_vector[1] +
                                spin_vector[2] * spin_vector[2]);

        spin_vectors_[atom_id][0] = spin_vector[0] / norm;
        spin_vectors_[atom_id][1] = spin_vector[1] / norm;
        spin_vectors_[atom_id][2] = spin_vector[2] / norm;
    }

    // Returns normalized spin vector [x, y, z] with ||spin|| = 1.0
    inline std::array<double, 3> generate_random_spin() const {
        double x, y, z, norm;
        do {
            x = Random::uniform_real(-1.0, 1.0);
            y = Random::uniform_real(-1.0, 1.0);
            z = Random::uniform_real(-1.0, 1.0);
            norm = x * x + y * y + z * z;
        } while (norm == 0.0 || norm > 1.0);
        norm = std::sqrt(norm);

        return {x / norm, y / norm, z / norm};
    }

    //
    inline int32_t select_random_magnetic_atom() const {
        return magnetic_atoms_[Random::uniform_int<int32_t>(
            0, magnetic_atoms_.size() - 1)];
    }

    // === INITIALIZING CONFIGURATIONS ===

    // Initializes all magnetic spins with random orientations
    void initialize_random();

    // Initializes all magnetic spins in ferromagnetic alignment (+z direction)
    void initialize_ferromagnetic();

    // Initializes magnetic spins in antiferromagnetic pattern (+-z alternating
    // direction)
    void initialize_antiferromagnetic();

    // Returns the configuration object
    inline const Config &get_config() const { return config_; }

    // Returns the geometry object
    inline const Geometry &get_geometry() const { return geometry_; }

  private:
    const Config &config_;
    const Geometry geometry_;

    std::vector<std::array<double, 3>> spin_vectors_;
    std::vector<bool> magnetic_mask_;

    mutable int32_t magnetic_count_ = 0;
    mutable std::vector<int32_t> magnetic_atoms_;
    mutable std::vector<int32_t> defect_atoms_;
    mutable bool cache_valid_ = false;

    void initialize_spins();
    void update_cache_data() const;
};

} // namespace lattice

#endif // ATOMS_HPP