#ifndef ATOMS_HPP
#define ATOMS_HPP

#include "../utils/random.hpp"
#include "geometry.hpp"

namespace lattice {

class Atoms {
  public:
    explicit Atoms(const Config &config);

    // === DEFECT MANAGEMENT ===

    // Sets the magnetic state of an atom (true = magnetic, false =
    // non-magnetic)
    inline void set_magnetic_state(int32_t atom_id, bool is_magnetic) {
        if (atom_id < 0 || atom_id >= geometry_.get_atom_count())
            throw std::out_of_range("Atom ID out of range");

        magnetic_mask_[atom_id] = is_magnetic;
    }

    // Returns the magnetic state of the specified atom
    inline bool get_magnetic_state(int32_t atom_id) const {
        if (atom_id < 0 || atom_id >= geometry_.get_atom_count())
            throw std::out_of_range("Atom ID out of range");

        return magnetic_mask_[atom_id];
    }

    // Returns a vector representing the magnetic mask for all atoms
    const std::vector<bool> get_magnetic_mask() { return magnetic_mask_; }

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
    inline std::array<double, 3> get_spin(int32_t atom_id) const {
        if (atom_id < 0 || atom_id >= geometry_.get_atom_count())
            throw std::out_of_range("Atom ID out of range");
        if (!magnetic_mask_[atom_id])
            throw std::invalid_argument("Atom is not magnetic");

        return spin_vectors_[atom_id];
    }

    // Sets the spin vector for the specified atom
    inline void set_spin(int32_t atom_id,
                         const std::array<double, 3> &spin_vector) {
        if (atom_id < 0 || atom_id >= geometry_.get_atom_count())
            throw std::out_of_range("Atom ID out of range");
        if (!magnetic_mask_[atom_id])
            throw std::invalid_argument(
                "Cannot set spin for non-magnetic atom");

        double norm = std::sqrt(spin_vector[0] * spin_vector[0] +
                                spin_vector[1] * spin_vector[1] +
                                spin_vector[2] * spin_vector[2]);

        spin_vectors_[atom_id] = {spin_vector[0] / norm, spin_vector[1] / norm,
                                  spin_vector[2] / norm};
    }

    // Returns normalized spin vector [x, y, z] with ||spin|| = 1.0
    inline std::array<double, 3> generate_random_spin() const {
        double x = Random::uniform_real(-1.0, 1.0);
        double y = Random::uniform_real(-1.0, 1.0);
        double z = Random::uniform_real(-1.0, 1.0);
        double norm = std::sqrt(x * x + y * y + z * z);

        return {x / norm, y / norm, z / norm};
    }

    //
    int32_t select_random_magnetic_atom() const;

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

    void initialize_spins();
};

} // namespace lattice

#endif // ATOMS_HPP