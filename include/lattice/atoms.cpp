#include "atoms.h"

lattice::Atoms::Atoms(const Config &config)
    : config_(config), geometry_(config)
{
    initialize_spins();
}

void lattice::Atoms::initialize_spins()
{
    int32_t num_atoms = geometry_.get_atom_count();
    spin_vectors_.resize(num_atoms);
    magnetic_mask_.resize(num_atoms, true);

    for (auto &spin : spin_vectors_)
        spin = {0.0, 0.0, 0.0};
}

// === DEFECT MANAGEMENT ===

void lattice::Atoms::set_random_defects(double defect_concentration)
{
    if (defect_concentration < 0 || defect_concentration > 1)
        throw std::invalid_argument("Defect concentration must be in [0,1]");

    int32_t num_atoms = geometry_.get_atom_count();
    std::fill(magnetic_mask_.begin(), magnetic_mask_.end(), true);

    int32_t defect_count = static_cast<int32_t>(defect_concentration * num_atoms);

    std::vector<int32_t> indices(num_atoms);
    for (int32_t i = 0; i < num_atoms; ++i)
        indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), Random::get_rng());

    for (int32_t i = 0; i < defect_count; ++i)
        magnetic_mask_[indices[i]] = false;
}

// === SYSTEM STATISTICS ===

int32_t lattice::Atoms::get_magnetic_count() const
{
    int32_t count = 0;
    for (bool is_magnetic : magnetic_mask_)
        count += is_magnetic;

    return count;
}

int32_t lattice::Atoms::get_defect_count() const
{
    int32_t count = 0;
    for (bool is_magnetic : magnetic_mask_)
        count += !is_magnetic;

    return count;
}

std::vector<int32_t> lattice::Atoms::get_magnetic_atoms() const
{
    std::vector<int32_t> magnetic_atoms;
    magnetic_atoms.reserve(magnetic_mask_.size());

    for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
        if (magnetic_mask_[i])
            magnetic_atoms.push_back(i);

    return magnetic_atoms;
}

std::vector<int32_t> lattice::Atoms::get_defect_atoms() const
{
    std::vector<int32_t> defect_atoms;
    defect_atoms.reserve(magnetic_mask_.size());

    for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
        if (!magnetic_mask_[i])
            defect_atoms.push_back(i);

    return defect_atoms;
}

// === WORKING WITH SPINS ===

int32_t lattice::Atoms::select_random_magnetic_atom() const
{
    const auto &magnetic_atoms = get_magnetic_atoms();
    return magnetic_atoms[Random::uniform_int<int32_t>(0, magnetic_atoms.size() - 1)];
}

// === INITIALIZING CONFIGURATIONS ===

void lattice::Atoms::initialize_random()
{
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i)
    {
        spin_vectors_[i] = generate_random_spin();
    }
}

void lattice::Atoms::initialize_ferromagnetic()
{
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i)
        spin_vectors_[i] = {0.0, 0.0, 1.0};
}

void lattice::Atoms::initialize_antiferromagnetic()
{
    for (int32_t i = 0; i < geometry_.get_atom_count(); ++i)
    {
        auto [cell_i, cell_j, atom_in_cell_id] = geometry_.get_cell_index(i);
        bool spin_up = (cell_i + cell_j + atom_in_cell_id) % 2 == 0;
        spin_vectors_[i] = {0.0, 0.0, spin_up ? 1.0 : -1.0};
    }
}