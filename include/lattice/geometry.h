#pragma once

#include <array>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include "../core/config.h"
#include "types.h"

namespace lattice
{

    class Geometry
    {
    public:
        explicit Geometry(const Config &config);

        int32_t get_linear_size() const { return linear_size_; }
        int32_t get_num_shells() const { return num_shells_; }
        int32_t get_num_positions() const { return num_positions_; }
        int32_t get_num_atoms() const { return num_atoms_; }
        double get_norm_a() const { return norm_a_; }
        double get_norm_b() const { return norm_b_; }
        CrystalType get_crystal_type() const { return crystal_type_; }
        BoundaryType get_boundary_type() const { return boundary_conditions_; }

        const std::vector<std::array<int32_t, 3>> &get_indexes() const { return indexes_; }
        const std::vector<std::array<double, 2>> &get_coordinates() const { return coordinates_; }
        const std::vector<std::vector<std::vector<int32_t>>> &get_neighbors() const { return neighbors_; }

        std::array<int32_t, 3> expand_idx(const int32_t idx) const;
        int32_t collapse_idx(const int32_t i, const int32_t j, const int32_t pos) const;
        std::array<double, 2> calculate_coordinate(int32_t idx) const;
        double calculate_distance(int32_t first_idx, int32_t second_idx) const;

    private:
        const int32_t linear_size_, num_shells_;
        int32_t num_positions_, num_atoms_;
        const double norm_a_, norm_b_;

        const CrystalType crystal_type_;
        const BoundaryType boundary_conditions_;

        std::vector<std::array<double, 2>> basis_;
        std::vector<std::array<double, 2>> translation_;
        std::vector<std::array<int32_t, 3>> indexes_;
        std::vector<std::array<double, 2>> coordinates_;
        std::vector<std::vector<std::vector<int32_t>>> neighbors_;

        static CrystalType parse_crystal_type(const std::string &crystal);
        static BoundaryType parse_boundary_type(const std::string &boundary);

        void initialize();
        void validate_parameters() const;
        void precomp_indexes();
        void precomp_coordinates();
        void precomp_neighbors();
    };

} // namespace lattice