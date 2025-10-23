#include "atoms.h"

lattice::Atoms::Atoms(const Config &config)
    : config_(config)
{
    initialize_geometry();
    initialize_spins();
}

void lattice::Atoms::initialize_geometry()
{
    geometry_ = std::make_unique<Geometry>(config_);
}

void lattice::Atoms::initialize_spins()
{
    int32_t num_atoms = geometry_->get_num_atoms();
    spins_.resize(num_atoms);
    magnetic_mask_.resize(num_atoms, true);

    for (auto &spin : spins_)
        spin = {0.0, 0.0, 0.0};
}

// === DEFECT MANAGEMENT ===

void lattice::Atoms::set_magnetic(int32_t atom_index, bool magnetic)
{
    if (atom_index < 0 || atom_index >= get_num_atoms())
        throw std::invalid_argument("Index is out of range num_atoms");

    magnetic_mask_[atom_index] = magnetic;
}

bool lattice::Atoms::is_magnetic(int32_t atom_index) const
{
    if (atom_index < 0 || atom_index >= get_num_atoms())
        throw std::invalid_argument("Index is out of range num_atoms");

    return magnetic_mask_[atom_index];
}

const std::vector<bool> lattice::Atoms::get_magnetic_mask() { return magnetic_mask_; }

void lattice::Atoms::set_random_defects(double concentration)
{
    if (concentration < 0 || concentration > 1)
        throw std::invalid_argument("Incorrect concentration value. Concentration must be in [0,1]");

    std::fill(magnetic_mask_.begin(), magnetic_mask_.end(), true);

    int32_t defect_count = static_cast<int32_t>(concentration * get_num_atoms());

    std::vector<int32_t> indices(get_num_atoms());
    for (int32_t i = 0; i < get_num_atoms(); ++i)
        indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), Random::get_rng());

    for (int32_t i = 0; i < defect_count; ++i)
        magnetic_mask_[indices[i]] = false;
}

// === SYSTEM STATISTICS ===

int32_t lattice::Atoms::count_magnetic_atoms() const
{
    int32_t count = 0;
    for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
        count += magnetic_mask_[i];

    return count;
}

int32_t lattice::Atoms::count_defects() const
{
    int32_t count = 0;
    for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
        count += !magnetic_mask_[i];

    return count;
}

std::vector<int32_t> lattice::Atoms::get_magnetic_atom_indices() const
{
    std::vector<int32_t> indices;
    indices.reserve(magnetic_mask_.size());

    for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
        if (magnetic_mask_[i])
            indices.push_back(i);

    return indices;
}

std::vector<int32_t> lattice::Atoms::get_defect_indices() const
{
    std::vector<int32_t> indices;
    indices.reserve(magnetic_mask_.size());

    for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
        if (!magnetic_mask_[i])
            indices.push_back(i);

    return indices;
}

// === WORKING WITH SPINS ===

std::array<double, 3> lattice::Atoms::get_spin(int32_t atom_index) const
{
    if (atom_index < 0 || atom_index >= get_num_atoms())
        throw std::invalid_argument("Atom index out of range");
    if (!magnetic_mask_[atom_index])
        throw std::invalid_argument("Atom is not magnetic");

    return spins_[atom_index];
}

void lattice::Atoms::set_spin(int32_t atom_index, const std::array<double, 3> &value)
{
    if (atom_index < 0 || atom_index >= get_num_atoms())
        throw std::invalid_argument("Atom index out of range");
    if (!magnetic_mask_[atom_index])
        throw std::invalid_argument("Cannot set spin for non-magnetic atom");

    double norm = std::sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
    spins_[atom_index] = {value[0] / norm, value[1] / norm, value[2] / norm};
}

void lattice::Atoms::flip_spins(const std::vector<int32_t> &indices)
{
    for (int32_t idx : indices)
    {
        if (idx < 0 || idx >= get_num_atoms())
            throw std::invalid_argument("Atom index out of range");
        if (!magnetic_mask_[idx])
            throw std::invalid_argument("Cannot flip spin of non-magnetic atom");

        spins_[idx] = {-spins_[idx][0], -spins_[idx][1], -spins_[idx][2]};
    }
}

// === INITIALIZING CONFIGURATIONS ===

void lattice::Atoms::random_initialize()
{
    for (int32_t i = 0; i < get_num_atoms(); ++i)
    {
        double x = Random::uniform_real(-1.0, 1.0);
        double y = Random::uniform_real(-1.0, 1.0);
        double z = Random::uniform_real(-1.0, 1.0);
        double norm = std::sqrt(x * x + y * y + z * z);

        spins_[i] = {x / norm, y / norm, z / norm};
    }
}

void lattice::Atoms::ferromagnetic_initialize()
{
    for (int32_t i = 0; i < get_num_atoms(); ++i)
        spins_[i] = {0.0, 0.0, 1.0};
}

void lattice::Atoms::antiferromagnetic_initialize()
{
    for (int32_t i = 0; i < get_num_atoms(); ++i)
        spins_[i] = {0.0, 0.0, (i % 2 == 0) ? 1.0 : -1.0};
}