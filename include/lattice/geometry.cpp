#include "geometry.hpp"
#include <algorithm>
#include <cmath>

lattice::Geometry::Geometry(const Config &config)
    : system_size_(config.get<int32_t>("lattice.system_size")),
      shell_count_(
          config.get<std::vector<double>>("physical.exchange_constants")
              .size()),
      lattice_constant_a_(config.get<double>("lattice.lattice_constant_a")),
      lattice_constant_b_(config.get<double>("lattice.lattice_constant_b")),
      crystal_type_(
          parse_crystal_type(config.get<std::string>("lattice.crystal_type"))),
      boundary_type_(parse_boundary_type(
          config.get<std::string>("lattice.boundary_type"))) {
    initialize_lattice();
    validate_parameters();
    compute_indices();
    compute_positions();
    compute_neighbors();
}

lattice::CrystalType
lattice::Geometry::parse_crystal_type(const std::string &crystal) {
    static const std::unordered_map<std::string, CrystalType> mapping = {
        {"rectangular", CrystalType::RECTANGULAR},
        {"c_rectangular", CrystalType::C_RECTANGULAR},
        {"triangular", CrystalType::TRIANGULAR},
        {"honeycomb", CrystalType::HONEYCOMB},
        {"kagome", CrystalType::KAGOME},
        {"lieb", CrystalType::LIEB}};

    if (auto it = mapping.find(crystal); it != mapping.end())
        return it->second;

    throw std::invalid_argument("Unknown crystal type: " + crystal);
}

lattice::BoundaryType
lattice::Geometry::parse_boundary_type(const std::string &boundary) {
    static const std::unordered_map<std::string, BoundaryType> mapping = {
        {"periodic", BoundaryType::PERIODIC}, {"hard", BoundaryType::HARD}};

    if (auto it = mapping.find(boundary); it != mapping.end())
        return it->second;
    throw std::invalid_argument("Unknown boundary type: " + boundary);
}

double lattice::Geometry::get_distance(int32_t first_atom_id,
                                       int32_t second_atom_id) const {
    if (first_atom_id < 0 || first_atom_id >= atom_count_ ||
        second_atom_id < 0 || second_atom_id > atom_count_)
        throw std::out_of_range(
            "Atom index out of range in calculate_distance");

    if (first_atom_id == second_atom_id)
        return 0.0;

    const double delta_x =
        atom_positions_[second_atom_id][0] - atom_positions_[first_atom_id][0];
    const double delta_y =
        atom_positions_[second_atom_id][1] - atom_positions_[first_atom_id][1];

    if (boundary_type_ == BoundaryType::PERIODIC) {
        double min_distance_sq = std::numeric_limits<double>::max();
        for (int32_t n = -1; n <= 1; ++n) {
            for (int32_t m = -1; m <= 1; ++m) {
                const double taco_dx =
                    (n * lattice_vectors_[0][0] + m * lattice_vectors_[1][0]) *
                        system_size_ +
                    delta_x;
                const double taco_dy =
                    (n * lattice_vectors_[0][1] + m * lattice_vectors_[1][1]) *
                        system_size_ +
                    delta_y;
                const double distance_sq =
                    taco_dx * taco_dx + taco_dy * taco_dy;
                if (distance_sq < min_distance_sq)
                    min_distance_sq = distance_sq;
            }
        }
        return std::sqrt(min_distance_sq);
    }

    return std::sqrt(delta_x * delta_x + delta_y * delta_y);
}

void lattice::Geometry::initialize_lattice() {
    const double SQRT3_2 = std::sqrt(3.0) / 2.0;

    switch (crystal_type_) {
    case CrystalType::RECTANGULAR:
        basis_count_ = 1;
        basis_vectors_ = {{0.0, 0.0}};
        lattice_vectors_ = {{lattice_constant_a_, 0.0},
                            {0.0, lattice_constant_b_}};
        break;

    case CrystalType::C_RECTANGULAR:
        basis_count_ = 2;
        basis_vectors_ = {
            {0.0, 0.0}, {0.5 * lattice_constant_a_, 0.5 * lattice_constant_b_}};
        lattice_vectors_ = {{lattice_constant_a_, 0.0},
                            {0.0, lattice_constant_b_}};
        break;

    case CrystalType::TRIANGULAR:
        basis_count_ = 1;
        basis_vectors_ = {{0.0, 0.0}};
        lattice_vectors_ = {
            {lattice_constant_a_, 0.0},
            {0.5 * lattice_constant_a_, SQRT3_2 * lattice_constant_b_}};
        break;

    case CrystalType::HONEYCOMB:
        basis_count_ = 2;
        basis_vectors_ = {
            {0.0, 0.0},
            {0.5 * lattice_constant_a_, SQRT3_2 / 3.0 * lattice_constant_b_}};
        lattice_vectors_ = {
            {lattice_constant_a_, 0.0},
            {0.5 * lattice_constant_a_, SQRT3_2 * lattice_constant_b_}};
        break;

    case CrystalType::KAGOME:
        basis_count_ = 3;
        basis_vectors_ = {
            {0.0, 0.0},
            {0.5 * lattice_constant_a_, 0.0},
            {0.25 * lattice_constant_a_, 0.5 * SQRT3_2 * lattice_constant_b_}};
        lattice_vectors_ = {
            {lattice_constant_a_, 0.0},
            {0.5 * lattice_constant_a_, SQRT3_2 * lattice_constant_b_}};
        break;

    case CrystalType::LIEB:
        basis_count_ = 3;
        basis_vectors_ = {{0.0, 0.0},
                          {0.5 * lattice_constant_a_, 0.0},
                          {0.0, 0.5 * lattice_constant_b_}};
        lattice_vectors_ = {{lattice_constant_a_, 0.0},
                            {0.0, lattice_constant_b_}};
        break;
    }

    atom_count_ = system_size_ * system_size_ * basis_count_;
}

void lattice::Geometry::validate_parameters() const {
    if (system_size_ <= 0)
        throw std::invalid_argument("Linear size must be positive, got: " +
                                    std::to_string(system_size_));

    if (shell_count_ <= 0)
        throw std::invalid_argument("Number of shells must be positive, got: " +
                                    std::to_string(shell_count_));

    if (lattice_constant_a_ <= 0 || lattice_constant_b_ <= 0)
        throw std::invalid_argument(
            "Lattice constants must be positive, got norm_a: " +
            std::to_string(lattice_constant_a_) +
            ", norm_b: " + std::to_string(lattice_constant_b_));

    if (lattice_vectors_.empty() || lattice_vectors_.size() != 2)
        throw std::runtime_error(
            "Translation vectors not properly initialized");

    const double total_system_x =
        (std::abs(lattice_vectors_[0][0]) + std::abs(lattice_vectors_[1][0])) *
        system_size_;
    const double total_system_y =
        (std::abs(lattice_vectors_[0][1]) + std::abs(lattice_vectors_[1][1])) *
        system_size_;

    const double max_system_distance =
        (boundary_type_ == BoundaryType::PERIODIC)
            ? std::sqrt(total_system_x * total_system_x +
                        total_system_y * total_system_y) /
                  2.0
            : std::sqrt(total_system_x * total_system_x +
                        total_system_y * total_system_y);

    if (static_cast<double>(shell_count_) > max_system_distance)
        throw std::invalid_argument(
            "Search radius (" + std::to_string(shell_count_) +
            ") exceeds system size (" + std::to_string(max_system_distance) +
            "). shells: " + std::to_string(shell_count_) +
            ", min_distance: " + std::to_string(shell_count_) +
            ", linear_size: " + std::to_string(system_size_));
}

void lattice::Geometry::compute_indices() {
    atom_indices_.resize(atom_count_);
    for (int32_t i = 0; i < atom_count_; ++i)
        atom_indices_[i] = get_cell_index(i);
}

void lattice::Geometry::compute_positions() {
    atom_positions_.resize(atom_count_);
    for (int32_t i = 0; i < atom_count_; ++i)
        atom_positions_[i] = get_atom_position(i);
}

void lattice::Geometry::compute_neighbors() {
    neighbor_table_.resize(atom_count_);
    for (auto &atom_neighbors : neighbor_table_)
        atom_neighbors.resize(shell_count_);

    std::vector<std::pair<int32_t, int32_t>> displacements;
    for (int32_t di = -shell_count_; di <= shell_count_; ++di)
        for (int32_t dj = -shell_count_; dj <= shell_count_; ++dj)
            if (di * di + dj * dj <= shell_count_ * shell_count_)
                displacements.emplace_back(di, dj);

    for (int32_t atom_id = 0; atom_id < atom_count_; ++atom_id) {
        std::vector<std::pair<int32_t, double>> atom_neighbors;

        auto [cell_i, cell_j, atom_in_cell_id] = atom_indices_[atom_id];

        for (const auto &[di, dj] : displacements) {
            for (int32_t new_atom_in_cell_id = 0;
                 new_atom_in_cell_id < basis_count_; ++new_atom_in_cell_id) {
                int32_t neighbor_id =
                    get_atom_id(cell_i + di, cell_j + dj, new_atom_in_cell_id);

                if (neighbor_id < 0 || neighbor_id >= atom_count_ ||
                    neighbor_id == atom_id)
                    continue;

                if (boundary_type_ == BoundaryType::HARD) {
                    auto [ni, nj, npos] = atom_indices_[neighbor_id];
                    if (ni < 0 || ni >= system_size_ || nj < 0 ||
                        nj >= system_size_)
                        continue;
                }

                double distance = get_distance(atom_id, neighbor_id);
                atom_neighbors.emplace_back(neighbor_id, distance);
            }
        }

        std::sort(
            atom_neighbors.begin(), atom_neighbors.end(),
            [](const auto &a, const auto &b) { return a.second < b.second; });

        if (atom_neighbors.empty())
            continue;

        double current_distance = atom_neighbors[0].second;
        int32_t current_shell = 0;

        for (const auto &[neighbor_id, distance] : atom_neighbors) {
            if (std::abs(distance - current_distance) >
                1e-5 * current_distance) {
                current_distance = distance;
                current_shell++;
                if (current_shell >= shell_count_)
                    break;
            }
            if (current_shell < shell_count_)
                neighbor_table_[atom_id][current_shell].push_back(neighbor_id);
        }
    }
}