#include "../../include/lattice/geometry.h"

Lattice::CrystalType Lattice::parse_crystal_type(const std::string &crystal)
{
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

Lattice::BoundaryType Lattice::parse_boundary_type(const std::string &boundary)
{
    static const std::unordered_map<std::string, BoundaryType> mapping = {
        {"periodic", BoundaryType::PERIODIC},
        {"hard", BoundaryType::HARD}};

    if (auto it = mapping.find(boundary); it != mapping.end())
        return it->second;
    throw std::invalid_argument("Unknown boundary type: " + boundary);
}

std::array<int32_t, 3> Lattice::expand_idx(const int32_t idx) const
{
    int32_t pos = idx % num_positions_;
    int32_t q = idx / num_positions_;
    int32_t i = q % linear_size_;
    int32_t j = q / linear_size_;

    return {i, j, pos};
}

int32_t Lattice::collapse_idx(const int32_t i, const int32_t j, const int32_t pos) const
{
    int32_t taco_i = (i % linear_size_ + linear_size_) % linear_size_;
    int32_t taco_j = (j % linear_size_ + linear_size_) % linear_size_;

    return (taco_i + taco_j * linear_size_) * num_positions_ + pos;
}

std::array<double, 2> Lattice::calculate_coordinate(int32_t idx) const
{
    auto [i, j, pos] = indexes_[idx];
    double x = i * translation_[0][0] + j * translation_[1][0] + basis_[pos][0];
    double y = i * translation_[0][1] + j * translation_[1][1] + basis_[pos][1];

    return {x, y};
}

double Lattice::calculate_distance(int32_t first_idx, int32_t second_idx) const
{
    if (first_idx < 0 || first_idx >= num_atoms_ || second_idx < 0 || second_idx >= num_atoms_)
        throw std::out_of_range("Atom index out of range in calculate_distance");

    if (first_idx == second_idx)
        return 0.0;

    const double delta_x = coordinates_[second_idx][0] - coordinates_[first_idx][0];
    const double delta_y = coordinates_[second_idx][1] - coordinates_[first_idx][1];

    if (boundary_conditions_ == BoundaryType::PERIODIC)
    {
        double min_distance_sq = std::numeric_limits<double>::max();
        for (int32_t n = -1; n <= 1; ++n)
        {
            for (int32_t m = -1; m <= 1; ++m)
            {
                const double taco_dx = n * translation_ax_ + m * translation_ay_ + delta_x;
                const double taco_dy = n * translation_bx_ + m * translation_by_ + delta_y;
                const double distance_sq = taco_dx * taco_dx + taco_dy * taco_dy;

                if (distance_sq < min_distance_sq)
                    min_distance_sq = distance_sq;
            }
        }
        return std::sqrt(min_distance_sq);
    }
    else
    {
        return std::sqrt(delta_x * delta_x + delta_y * delta_y);
    }
}

void Lattice::initialize()
{
    const double SQRT3_2 = std::sqrt(3.0) / 2.0;

    switch (crystal_type_)
    {
    case CrystalType::RECTANGULAR:
        num_positions_ = 1;
        basis_ = {{0.0, 0.0}};
        translation_ = {
            {norm_a_, 0.0},
            {0.0, norm_b_}};
        break;

    case CrystalType::C_RECTANGULAR:
        num_positions_ = 2;
        basis_ = {
            {0.0, 0.0},
            {0.5 * norm_a_, 0.5 * norm_b_}};
        translation_ = {
            {norm_a_, 0.0},
            {0.0, norm_b_}};
        break;

    case CrystalType::TRIANGULAR:
        num_positions_ = 1;
        basis_ = {{0.0, 0.0}};
        translation_ = {
            {norm_a_, 0.0},
            {0.5 * norm_a_, SQRT3_2 * norm_b_}};
        break;

    case CrystalType::HONEYCOMB:
        num_positions_ = 2;
        basis_ = {{0.0, 0.0},
                  {0.5 * norm_a_, SQRT3_2 / 3.0 * norm_b_}};
        translation_ = {
            {norm_a_, 0.0},
            {0.5 * norm_a_, SQRT3_2 * norm_b_}};
        break;

    case CrystalType::KAGOME:
        num_positions_ = 3;
        basis_ = {
            {0.0, 0.0},
            {0.5 * norm_a_, 0.0},
            {0.25 * norm_a_, 0.5 * SQRT3_2 * norm_b_}};
        translation_ = {
            {norm_a_, 0.0},
            {0.5 * norm_a_, SQRT3_2 * norm_b_}};

        break;

    case CrystalType::LIEB:
        num_positions_ = 3;
        basis_ = {
            {0.0, 0.0},
            {0.5 * norm_a_, 0.0},
            {0.0, 0.5 * norm_b_}};
        translation_ = {
            {norm_a_, 0.0},
            {0.0, norm_b_}};
        break;
    }
    num_atoms_ = linear_size_ * linear_size_ * num_positions_;

    translation_ax_ = translation_[0][0] * linear_size_;
    translation_ay_ = translation_[0][1] * linear_size_;
    translation_bx_ = translation_[1][0] * linear_size_;
    translation_by_ = translation_[1][1] * linear_size_;
}

void Lattice::validate_parameters() const
{
    const double lattice_diagonal = std::sqrt(translation_ax_ * translation_ax_ + translation_by_ * translation_by_);
    const double max_possible_distance = (boundary_conditions_ == BoundaryType::PERIODIC)
                                             ? lattice_diagonal / 2.0
                                             : lattice_diagonal;

    if (static_cast<double>(num_shells_) > max_possible_distance)
        throw std::invalid_argument("Number of shells (" +
                                    std::to_string(num_shells_) +
                                    ") exceeds maximum possible distance (" +
                                    std::to_string(max_possible_distance) +
                                    ") in the lattice");
}

void Lattice::precomp_indexes()
{
    indexes_.resize(num_atoms_);
    for (int32_t i = 0; i < num_atoms_; ++i)
        indexes_[i] = expand_idx(i);
}

void Lattice::precomp_coordinates()
{
    coordinates_.resize(num_atoms_);
    for (int32_t i = 0; i < num_atoms_; ++i)
        coordinates_[i] = calculate_coordinate(i);
}

void Lattice::precomp_neighbors()
{
    neighbors_.resize(num_atoms_);
    for (int32_t i = 0; i < num_atoms_; ++i)
        neighbors_[i].resize(num_shells_);

    std::vector<bool> atom_visited(num_atoms_, false);

    for (int32_t atom_idx = 0; atom_idx < num_atoms_; ++atom_idx)
    {
        std::fill(atom_visited.begin(), atom_visited.end(), false);
        atom_visited[atom_idx] = true;

        std::vector<std::pair<int32_t, double>> atom_neighbors;

        for (int32_t di = -num_shells_; di <= num_shells_; ++di)
        {
            for (int32_t dj = -num_shells_; dj <= num_shells_; ++dj)
            {
                if (di * di + dj * dj > num_shells_ * num_shells_)
                    continue;

                auto [i, j, pos] = indexes_[atom_idx];

                for (int32_t new_pos = 0; new_pos < num_positions_; ++new_pos)
                {
                    int32_t neighbor_idx = collapse_idx(i + di, j + dj, new_pos);

                    if (neighbor_idx < 0 || neighbor_idx >= num_atoms_)
                        continue;
                    if (atom_visited[neighbor_idx])
                        continue;
                    if (neighbor_idx == atom_idx)
                        continue;
                    if (boundary_conditions_ == BoundaryType::HARD)
                    {
                        auto [ni, nj, npos] = indexes_[neighbor_idx];
                        if (ni < 0 || ni >= linear_size_ || nj < 0 || nj >= linear_size_)
                            continue;
                    }

                    atom_visited[neighbor_idx] = true;
                    double distance = calculate_distance(atom_idx, neighbor_idx);
                    atom_neighbors.emplace_back(neighbor_idx, distance);
                }
            }
        }

        std::sort(atom_neighbors.begin(), atom_neighbors.end(), [](const auto &a, const auto &b)
                  { return a.second < b.second; });

        if (atom_neighbors.empty())
            continue;

        double current_distance = atom_neighbors[0].second;
        int32_t current_shell = 0;

        for (const auto &[neighbor_idx, distance] : atom_neighbors)
        {
            if (std::abs(distance - current_distance) > 1e-5)
            {
                current_distance = distance;
                current_shell++;
                if (current_shell >= num_shells_)
                    break;
            }
            if (current_shell < num_shells_)
                neighbors_[atom_idx][current_shell].push_back(neighbor_idx);
        }
    }
}