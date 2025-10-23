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

        // Returns the linear size of the system (number of unit cells along each dimension)
        int32_t get_system_size() const { return system_size_; }

        // Returns the number of coordination shells
        int32_t get_shell_count() const { return shell_count_; }

        // Returns the number of atomic positions in the unit cell basis
        int32_t get_basis_count() const { return basis_count_; }

        // Returns the total number of atoms in the system
        int32_t get_atom_count() const { return atom_count_; }

        // Returns the lattice constant along the a-direction
        double get_lattice_constant_a() const { return lattice_constant_a_; }

        // Returns the lattice constant along the b-direction
        double get_lattice_constant_b() const { return lattice_constant_b_; }

        // Returns the crystal structure type
        CrystalType get_crystal_type() const { return crystal_type_; }

        // Returns the boundary conditions type
        BoundaryType get_boundary_type() const { return boundary_type_; }

        // Returns all atom indices as [cell_i, cell_j, atom_in_cell_id] triplets
        const std::vector<std::array<int32_t, 3>> &get_atom_indices() const { return atom_indices_; }

        // Returns all atom coordinates as [x, y] pairs
        const std::vector<std::array<double, 2>> &get_atom_positions() const { return atom_positions_; }

        // Returns the neighbor table: [atom_id][shell][neighbor_index]
        const std::vector<std::vector<std::vector<int32_t>>> &get_neighbor_table() const { return neighbor_table_; }

        // Converts atom ID to cell indices [cell_i, cell_j, atom_in_cell_id]
        std::array<int32_t, 3> get_cell_index(const int32_t atom_id) const;

        // Converts cell indices to atom ID
        int32_t get_atom_id(const int32_t cell_i, const int32_t cell_j, const int32_t atom_in_cell_id) const;

        // Returns the Cartesian coordinates [x, y] of the specified atom
        std::array<double, 2> get_atom_position(int32_t atom_id) const;

        // Calculates the distance between two atoms
        double get_distance(int32_t first_atom_id, int32_t second_atom_id) const;

    private:
        const int32_t system_size_, shell_count_;
        int32_t basis_count_, atom_count_;
        const double lattice_constant_a_, lattice_constant_b_;

        const CrystalType crystal_type_;
        const BoundaryType boundary_type_;

        std::vector<std::array<double, 2>> basis_vectors_;
        std::vector<std::array<double, 2>> lattice_vectors_;
        std::vector<std::array<int32_t, 3>> atom_indices_;
        std::vector<std::array<double, 2>> atom_positions_;
        std::vector<std::vector<std::vector<int32_t>>> neighbor_table_;

        static CrystalType parse_crystal_type(const std::string &crystal);
        static BoundaryType parse_boundary_type(const std::string &boundary);

        void initialize_lattice();
        void validate_parameters() const;
        void compute_indices();
        void compute_positions();
        void compute_neighbors();
    };

} // namespace lattice