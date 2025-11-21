#include "atoms.hpp"
#include <algorithm>
#include <array>
#include <cstdint>

lattice::Atoms::Atoms(const Config &config)
    : config_(config), geometry_(config) {
    initialize_spins();
}

void lattice::Atoms::initialize_spins() {
    int32_t num_atoms = geometry_.get_atom_count();
    spin_vectors_.resize(num_atoms);
    magnetic_mask_.resize(num_atoms, true);
    magnetic_atoms_.reserve(num_atoms);
    defect_atoms_.reserve(num_atoms);

    update_cache_data();
}

// === DEFECT MANAGEMENT ===

void lattice::Atoms::set_random_defects(double defect_concentration) {
    if (defect_concentration < 0 || defect_concentration > 1)
        throw std::invalid_argument("Defect concentration must be in [0,1]");

    int32_t num_atoms = geometry_.get_atom_count();
    std::fill(magnetic_mask_.begin(), magnetic_mask_.end(), true);

    int32_t defect_count =
        static_cast<int32_t>(defect_concentration * num_atoms);

    // Fisher–Yates shuffle
    for (int32_t i = num_atoms - 1; i > 0; --i) {
        int32_t j = Random::uniform_int(0, i);

        bool temp = magnetic_mask_[i];
        magnetic_mask_[i] = magnetic_mask_[j];
        magnetic_mask_[j] = temp;
    }

    for (int32_t i = 0; i < defect_count; ++i)
        magnetic_mask_[i] = false;

    update_cache_data();
}

// === INITIALIZING CONFIGURATIONS ===

void lattice::Atoms::initialize_random() {
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i) {
        spin_vectors_[i] = generate_random_spin();
    }
}

void lattice::Atoms::initialize_ferromagnetic() {
    const std::array<double, 3> spin = {0.0, 0.0, 1.0};
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i)
        spin_vectors_[i] = spin;
}

void lattice::Atoms::initialize_antiferromagnetic() {
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i) {
        auto [cell_i, cell_j, atom_in_cell_id] = geometry_.get_cell_index(i);
        bool spin_up = (cell_i + cell_j + atom_in_cell_id) % 2 == 0;
        spin_vectors_[i] = {0.0, 0.0, spin_up ? 1.0 : -1.0};
    }
}

void lattice::Atoms::update_cache_data() const {
    int32_t num_atoms = geometry_.get_atom_count();

    magnetic_count_ = 0;
    magnetic_atoms_.clear();
    defect_atoms_.clear();

    magnetic_atoms_.reserve(num_atoms);
    defect_atoms_.reserve(num_atoms);

    for (int32_t i = 0; i < num_atoms; ++i) {
        if (magnetic_mask_[i]) {
            magnetic_atoms_.push_back(i);
            ++magnetic_count_;
        } else
            defect_atoms_.push_back(i);
    }

    magnetic_atoms_.shrink_to_fit();
    defect_atoms_.shrink_to_fit();
    cache_valid_ = true;
}