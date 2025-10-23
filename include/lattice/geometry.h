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
        inline int32_t get_system_size() const { return system_size_; }

        // Returns the number of coordination shells
        inline int32_t get_shell_count() const { return shell_count_; }

        // Returns the number of atomic positions in the unit cell basis
        inline int32_t get_basis_count() const { return basis_count_; }

        // Returns the total number of atoms in the system
        inline int32_t get_atom_count() const { return atom_count_; }

        // Returns the lattice constant along the a-direction
        inline double get_lattice_constant_a() const { return lattice_constant_a_; }

        // Returns the lattice constant along the b-direction
        inline double get_lattice_constant_b() const { return lattice_constant_b_; }

        // Returns the crystal structure type
        inline CrystalType get_crystal_type() const { return crystal_type_; }

        // Returns the boundary conditions type
        inline BoundaryType get_boundary_type() const { return boundary_type_; }

        // Returns all atom indices as [cell_i, cell_j, atom_in_cell_id] triplets
        inline const std::vector<std::array<int32_t, 3>> &get_atom_indices() const { return atom_indices_; }

        // Returns all atom coordinates as [x, y] pairs
        inline const std::vector<std::array<double, 2>> &get_atom_positions() const { return atom_positions_; }

        // Returns the neighbor table: [atom_id][shell][neighbor_index]
        inline const std::vector<std::vector<std::vector<int32_t>>> &get_neighbor_table() const { return neighbor_table_; }

        // Converts atom ID to cell indices [cell_i, cell_j, atom_in_cell_id]
        inline std::array<int32_t, 3> get_cell_index(const int32_t atom_id) const
        {
            int32_t atom_in_cell_id = atom_id % basis_count_;
            int32_t temp = atom_id / basis_count_;
            int32_t cell_i = temp % system_size_;
            int32_t cell_j = temp / system_size_;

            return {cell_i, cell_j, atom_in_cell_id};
        }

        // Converts cell indices to atom ID
        inline int32_t get_atom_id(const int32_t cell_i, const int32_t cell_j, const int32_t atom_in_cell_id) const
        {
            int32_t taco_i = (cell_i % system_size_ + system_size_) % system_size_;
            int32_t taco_j = (cell_j % system_size_ + system_size_) % system_size_;

            return (taco_i + taco_j * system_size_) * basis_count_ + atom_in_cell_id;
        }

        // Returns the Cartesian coordinates [x, y] of the specified atom
        inline std::array<double, 2> get_atom_position(int32_t atom_id) const
        {
            auto [cell_i, cell_j, atom_in_cell_id] = atom_indices_[atom_id];
            double x = cell_i * lattice_vectors_[0][0] + cell_j * lattice_vectors_[1][0] + basis_vectors_[atom_in_cell_id][0];
            double y = cell_i * lattice_vectors_[0][1] + cell_j * lattice_vectors_[1][1] + basis_vectors_[atom_in_cell_id][1];

            return {x, y};
        }

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