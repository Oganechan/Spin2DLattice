#pragma once

#include <array>
#include <vector>
#include "config.h"
#include "random.h"

class Lattice
{
public:
    enum class CrystalType
    {
        Rectangular,
        Square,
        Hexagonal,
    };
    enum class BoundaryType
    {
        Periodic,
        Hard
    };

    explicit Lattice(const Config &config)
        : linear_size_(config.get<size_t>("lattice.linear_size")),
          num_shells_(config.get<size_t>("lattice.num_shells")),
          num_positions_(config.get<size_t>("lattice.num_positions")),
          num_atoms_(linear_size_ * linear_size_ * num_positions_),
          norm_a_(config.get<double>("lattice.norm_a", 1.0)),
          norm_b(config.get<double>("lattice.norm_b", 1.0)),
          crystal_type_(parse_crystal_type(config.get<std::string>("lattice.crystal_type", "square"))),
          boundary_conditions_(parse_boundary_type(config.get<std::string>("lattice.boundary_conditions", "periodic")))
    {
        initialize();
        cached_neighbors_ = generate_neighbors();
    }

private:
    const size_t linear_size_, num_shells_, num_positions_, num_atoms_;
    const double norm_a_, norm_b;

    const CrystalType crystal_type_;
    static CrystalType parse_crystal_type(const std::string &type);

    const BoundaryType boundary_conditions_;
    static BoundaryType parse_boundary_type(const std::string &boundary);

    std::vector<std::array<double, 2>> basis_;
    std::vector<std::array<double, 2>> translation_;
    std::vector<std::array<size_t, 3>> indexes_;

    mutable std::vector<std::vector<std::vector<size_t>>> cached_neighbors_;
    std::vector<std::vector<std::vector<size_t>>> generate_neighbors() const;

    void initialize();
};