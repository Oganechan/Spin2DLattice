#include "../include/lattice.h"

Lattice::CrystalType Lattice::parse_crystal_type(const std::string &type)
{
    static const std::unordered_map<std::string, CrystalType> mapping = {
        {"rectangular", CrystalType::Rectangular},
        {"square", CrystalType::Square},
        {"hexagonal", CrystalType::Hexagonal}};

    if (auto it = mapping.find(type); it != mapping.end())
        return it->second;

    throw std::invalid_argument("Unknown crystal type: " + type);
}

Lattice::BoundaryType Lattice::parse_boundary_type(const std::string &boundary)
{
    static const std::unordered_map<std::string, BoundaryType> mapping = {
        {"periodic", BoundaryType::Periodic},
        {"hard", BoundaryType::Hard}};

    if (auto it = mapping.find(boundary); it != mapping.end())
        return it->second;
    throw std::invalid_argument("Unknown boundary type: " + boundary);
}

void Lattice::initialize()
{
}

std::vector<std::vector<std::vector<size_t>>> Lattice::generate_neighbors() const
{
}