#include "../include/lattice.h"

Lattice::CrystalType Lattice::parse_crystal_type(const std::string &crystal)
{
    static const std::unordered_map<std::string, CrystalType> mapping = {
        {"rectangular", CrystalType::RECTANGULAR},
        {"triangular", CrystalType::TRIANGULAR},
        {"honeycomb", CrystalType::HONEYCOMB},
        {"kagome", CrystalType::KAGOME},
        {"lieb", CrystalType::LIEB},
        {"checkerboard", CrystalType::CHECKERBOARD}};

    if (auto it = mapping.find(crystal); it != mapping.end())
        return it->second;

    throw std::invalid_argument("Unknown crystal type: " + crystal);
}

Lattice::BoundaryType Lattice::parse_boundary_type(const std::string &boundary)
{
    static const std::unordered_map<std::string, BoundaryType> mapping = {
        {"periodic", BoundaryType::PERIODIC},
        {"hard", BoundaryType::HARD}};

    if (auto it = mapping.find(boundary); it != mapping.end())
        return it->second;
    throw std::invalid_argument("Unknown boundary type: " + boundary);
}

void Lattice::initialize()
{
    switch (crystal_type_)
    {
    case CrystalType::RECTANGULAR:
        basis_ = {{0.0, 0.0}};
        translation_ = {
            {norm_a_, 0.0},
            {0.0, norm_b_}};
        num_positions_ = 1;
        break;

    case CrystalType::TRIANGULAR:
        basis_ = {{0.0, 0.0}};
        translation_ = {
            {norm_a_, 0.0},
            {0.5 * norm_a_, std::sqrt(3.0) / 2.0 * norm_b_}};
        num_positions_ = 1;
        break;

    case CrystalType::HONEYCOMB:
        basis_ = {{0.0, 0.0},
                  {1.0 / 3.0, 1.0 / 3.0}};
        translation_ = {
            {norm_a_, 0.0},
            {0.5 * norm_a_, std::sqrt(3.0) / 2.0 * norm_b_}};
        num_positions_ = 2;
        break;

    case CrystalType::KAGOME:
        basis_ = {
            {0.0, 0.0},
            {0.5, 0.0},
            {0.25, 0.5 * std::sqrt(3.0) / 2.0}};
        translation_ = {
            {norm_a_, 0.0},
            {0.5 * norm_a_, std::sqrt(3.0) / 2.0 * norm_b_}};
        num_positions_ = 3;
        break;

    case CrystalType::LIEB:
        basis_ = {
            {0.0, 0.0},
            {0.5, 0.0},
            {0.0, 0.5}};
        translation_ = {
            {norm_a_, 0.0},
            {0.0, norm_b_}};
        num_positions_ = 3;
        break;

    case CrystalType::CHECKERBOARD:
        basis_ = {
            {0.0, 0.0},
            {0.5, 0.5}};
        translation_ = {
            {norm_a_, 0.0},
            {0.0, norm_b_}};
        num_positions_ = 2;
        break;
    }

    num_atoms_ = linear_size_ * linear_size_ * num_positions_;
    indexes_.resize(num_atoms_);

    for (size_t i = 0; i < num_atoms_; ++i)
        indexes_[i] = expand_idx(i);
}

std::array<size_t, 3> Lattice::expand_idx(const size_t idx)
{
    const size_t pos = idx % num_positions_;
    const size_t q = idx / num_positions_;
    const size_t i = q % linear_size_;
    const size_t j = q / linear_size_;

    return {i, j, pos};
}

size_t Lattice::collapse_idx(const int i, const int j, const size_t pos)
{
    size_t taco_i = (i % static_cast<int>(linear_size_) + linear_size_) % linear_size_;
    size_t taco_j = (j % static_cast<int>(linear_size_) + linear_size_) % linear_size_;

    return (taco_i + taco_j * linear_size_) * num_positions_ + pos;
}

std::vector<std::vector<std::vector<size_t>>> Lattice::generate_neighbors() const
{
}