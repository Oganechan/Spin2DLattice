#pragma once

#include <array>
#include <vector>
#include <cmath>
#include "config.h"
#include "random.h"

class Lattice
{
public:
    enum class CrystalType
    {
        RECTANGULAR,
        TRIANGULAR,
        HONEYCOMB,
        KAGOME,
        LIEB,
        CHECKERBOARD
    };
    enum class BoundaryType
    {
        PERIODIC,
        HARD
    };

    explicit Lattice(const Config &config)
        : linear_size_(config.get<size_t>("lattice.linear_size")),
          num_shells_(config.get<size_t>("lattice.num_shells")),
          norm_a_(config.get<double>("lattice.norm_a", 1.0)),
          norm_b_(config.get<double>("lattice.norm_b", 1.0)),
          crystal_type_(parse_crystal_type(config.get<std::string>("lattice.crystal", "rectangular"))),
          boundary_conditions_(parse_boundary_type(config.get<std::string>("lattice.boundary", "periodic")))
    {
        initialize();
        cached_neighbors_ = generate_neighbors();
    }

private:
    const size_t linear_size_, num_shells_;
    size_t num_positions_, num_atoms_;
    const double norm_a_, norm_b_;

    const CrystalType crystal_type_;
    static CrystalType parse_crystal_type(const std::string &crystal);

    const BoundaryType boundary_conditions_;
    static BoundaryType parse_boundary_type(const std::string &boundary);

    std::vector<std::array<double, 2>> basis_;
    std::vector<std::array<double, 2>> translation_;
    std::vector<std::array<size_t, 3>> indexes_;

    mutable std::vector<std::vector<std::vector<size_t>>> cached_neighbors_;
    std::vector<std::vector<std::vector<size_t>>> generate_neighbors() const;

    void initialize();

    std::array<size_t, 3> expand_idx(const size_t idx);
    size_t collapse_idx(const int i, const int j, const size_t pos);
};