#include "atoms.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

lattice::Atoms::Atoms(const Config &config) : geometry_(config) {
    initialize_spins();
}

// Randomly sets atoms as non-magnetic atoms with given concentration
void lattice::Atoms::set_random_defects(double defect_concentration) {
    if (defect_concentration < 0 || defect_concentration > 1)
        throw std::invalid_argument("Defect concentration must be in [0,1]");

    const int32_t atom_count = geometry_.get_atom_count();
    std::fill(magnetic_mask_.begin(), magnetic_mask_.end(), true);

    const int32_t defect_count = int32_t(defect_concentration * atom_count);

    std::vector<int32_t> atom_indices(atom_count);
    std::iota(atom_indices.begin(), atom_indices.end(), 0);

    // Fisher–Yates shuffle
    for (int32_t i = atom_count - 1; i > 0; --i) {
        int32_t j = Random::uniform_int(0, i);
        std::swap(atom_indices[i], atom_indices[j]);
    }

    for (int32_t i = 0; i < defect_count; ++i)
        magnetic_mask_[atom_indices[i]] = false;

    update_cache_data();
}

std::array<double, 3> lattice::Atoms::small_rotate_spin(int32_t atom_id) {
    const auto &spin = get_spin(atom_id);

    static double max_angel = 30.0 * (M_PI / 180.0);
    double phi = std::atan2(spin[1], spin[0]) +
                 Random::uniform_real<double>(-max_angel, max_angel);
    double theta = std::acos(spin[2]) +
                   Random::uniform_real<double>(-max_angel, max_angel);

    double x = std::sin(theta) * std::cos(phi);
    double y = std::sin(theta) * std::sin(phi);
    double z = std::cos(theta);

    return {x, y, z};
}

// Initializes all magnetic spins with random orientations
void lattice::Atoms::initialize_random() {
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i) {
        spin_vectors_[i] = generate_random_spin();
    }
}

// Initializes all magnetic spins in ferromagnetic alignment (+z direction)
void lattice::Atoms::initialize_ferromagnetic() {
    const std::array<double, 3> spin = {0.0, 0.0, 1.0};
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i)
        spin_vectors_[i] = spin;
}

// Initializes magnetic spins in antiferromagnetic pattern (+-z alternating
// direction)
void lattice::Atoms::initialize_antiferromagnetic() {
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i) {
        auto [cell_i, cell_j, atom_in_cell_id] = geometry_.get_cell_index(i);
        bool spin_up = (cell_i + cell_j + atom_in_cell_id) % 2 == 0;
        spin_vectors_[i] = {0.0, 0.0, spin_up ? 1.0 : -1.0};
    }
}

void lattice::Atoms::initialize_spins() {
    int32_t atom_count = geometry_.get_atom_count();

    spin_vectors_.resize(atom_count);
    magnetic_mask_.resize(atom_count, true);
    magnetic_atoms_.reserve(atom_count);
    defect_atoms_.reserve(atom_count);

    update_cache_data();
}

void lattice::Atoms::update_cache_data() const {
    int32_t atom_count = geometry_.get_atom_count();

    magnetic_count_ = 0;
    magnetic_atoms_.clear();
    defect_atoms_.clear();

    magnetic_atoms_.reserve(atom_count);
    defect_atoms_.reserve(atom_count);

    for (int32_t i = 0; i < atom_count; ++i) {
        if (magnetic_mask_[i]) {
            magnetic_atoms_.push_back(i);
            ++magnetic_count_;
        } else
            defect_atoms_.push_back(i);
    }

    magnetic_atoms_.shrink_to_fit();
    defect_atoms_.shrink_to_fit();
}