#include "atoms.hpp"
#include "spin.hpp"
#include "types.hpp"
#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

lattice::Atoms::Atoms(const Config &config)
    : geometry_(config), spin_model_(parse_spin_model(
                             config.get<std::string>("physics.model_type"))) {
    initialize_spins();
}

lattice::SpinModel lattice::Atoms::parse_spin_model(const std::string &model) {
    static const std::unordered_map<std::string, SpinModel> mapping = {
        {"Ising", SpinModel::ISING},
        {"XY", SpinModel::XY},
        {"Heisenberg", SpinModel::HEISENBERG}};

    if (auto it = mapping.find(model); it != mapping.end())
        return it->second;

    throw std::invalid_argument("Unknowing spin model:" + model);
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

std::unique_ptr<lattice::BaseSpin>
lattice::Atoms::generate_random_spin() const {
    auto spin = spin_factory_.create(spin_model_);
    spin->randomize();
    return spin;
}

// Initializes all magnetic spins with random orientations
void lattice::Atoms::initialize_random() {
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i) {
        spins_[i]->randomize();
    }
}

void lattice::Atoms::initialize_spins() {
    int32_t atom_count = geometry_.get_atom_count();

    spins_.resize(atom_count);
    magnetic_mask_.resize(atom_count, true);
    magnetic_atoms_.reserve(atom_count);
    defect_atoms_.reserve(atom_count);

    for (int32_t i = 0; i < atom_count; ++i)
        spins_[i] = spin_factory_.create(spin_model_);

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