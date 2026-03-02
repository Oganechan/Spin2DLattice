#ifndef ATOMS_HPP
#define ATOMS_HPP

#include "../utils/random.hpp"
#include "geometry.hpp"
#include "spin.hpp"
#include "types.hpp"
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace lattice {

class Atoms {
  public:
    explicit Atoms(const Config &config);

    void set_random_defects(double defect_concentration);
    std::unique_ptr<BaseSpin> generate_random_spin() const;
    void initialize_random();

    // Sets the magnetic state of an atom (true = magnetic, false =
    // non-magnetic)
    inline void set_magnetic_state(int32_t atom_id, bool is_magnetic) {
        magnetic_mask_[atom_id] = is_magnetic;
        update_cache_data();
    }

    // Returns the magnetic state of the specified atom
    inline bool get_magnetic_state(int32_t atom_id) const {
        return magnetic_mask_[atom_id];
    }

    // Returns a vector representing the magnetic mask for all atoms
    inline const std::vector<bool> &get_magnetic_mask() const {
        return magnetic_mask_;
    }

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

    // Returns the spin vector of the specified atom
    inline std::array<double, 3> get_spin_components(int32_t atom_id) const {
        return spins_[atom_id]->get_components();
    }

    // Returns the spin object of the specified atom
    inline const BaseSpin &get_spin(int32_t atom_id) const {
        return *spins_[atom_id];
    }

    // Sets the spin vector for the specified atom
    inline void set_spin(int32_t atom_id, const BaseSpin &new_spin) {

        spins_[atom_id]->replace(new_spin);
    }

    // Returns random magnetic atom
    inline int32_t select_random_magnetic_atom() const {
        return magnetic_atoms_[rng_.uniform_int<int32_t>(
            0, magnetic_atoms_.size() - 1)];
    }

    // Returns the geometry object
    inline const Geometry &get_geometry() const { return geometry_; }

    // Returns spin type
    inline const SpinModel get_spin_model() const { return spin_model_; }

  private:
    const Geometry geometry_;
    const SpinModel spin_model_;
    Random &rng_ = thread_local_random();

    std::vector<std::unique_ptr<BaseSpin>> spins_;
    std::vector<bool> magnetic_mask_;

    mutable int32_t magnetic_count_;
    mutable std::vector<int32_t> magnetic_atoms_;
    mutable std::vector<int32_t> defect_atoms_;

    static SpinModel parse_spin_model(const std::string &model);

    void initialize_spins();
    void update_cache_data() const;
};

} // namespace lattice

#endif // ATOMS_HPP