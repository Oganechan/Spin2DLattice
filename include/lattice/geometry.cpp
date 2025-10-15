#include "geometry.h"

lattice::Geometry::Geometry(const Config &config)
    : linear_size_(config.get<int32_t>("lattice.linear_size")),
      num_shells_(config.get<std::vector<double>>("physical.exchange_interaction").size()),
      norm_a_(config.get<double>("lattice.norm_a")),
      norm_b_(config.get<double>("lattice.norm_b")),
      crystal_type_(parse_crystal_type(config.get<std::string>("lattice.crystal_type"))),
      boundary_conditions_(parse_boundary_type(config.get<std::string>("lattice.boundary_conditions")))
{
    initialize();
    validate_parameters();
    precomp_indexes();
    precomp_coordinates();
    precomp_neighbors();
}

lattice::CrystalType lattice::Geometry::parse_crystal_type(const std::string &crystal)
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

lattice::BoundaryType lattice::Geometry::parse_boundary_type(const std::string &boundary)
{
    static const std::unordered_map<std::string, BoundaryType> mapping = {
        {"periodic", BoundaryType::PERIODIC},
        {"hard", BoundaryType::HARD}};

    if (auto it = mapping.find(boundary); it != mapping.end())
        return it->second;
    throw std::invalid_argument("Unknown boundary type: " + boundary);
}
std::array<int32_t, 3> lattice::Geometry::expand_idx(const int32_t idx) const
{
    int32_t pos = idx % num_positions_;
    int32_t q = idx / num_positions_;
    int32_t i = q % linear_size_;
    int32_t j = q / linear_size_;

    return {i, j, pos};
}

int32_t lattice::Geometry::collapse_idx(const int32_t i, const int32_t j, const int32_t pos) const
{
    int32_t taco_i = (i % linear_size_ + linear_size_) % linear_size_;
    int32_t taco_j = (j % linear_size_ + linear_size_) % linear_size_;

    return (taco_i + taco_j * linear_size_) * num_positions_ + pos;
}

std::array<double, 2> lattice::Geometry::calculate_coordinate(int32_t idx) const
{
    auto [i, j, pos] = indexes_[idx];
    double x = i * translation_[0][0] + j * translation_[1][0] + basis_[pos][0];
    double y = i * translation_[0][1] + j * translation_[1][1] + basis_[pos][1];

    return {x, y};
}

double lattice::Geometry::calculate_distance(int32_t first_idx, int32_t second_idx) const
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
                const double taco_dx = (n * translation_[0][0] + m * translation_[1][0]) * linear_size_ + delta_x;
                const double taco_dy = (n * translation_[0][1] + m * translation_[1][1]) * linear_size_ + delta_y;
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

void lattice::Geometry::initialize()
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
}

void lattice::Geometry::validate_parameters() const
{
    if (linear_size_ <= 0)
        throw std::invalid_argument("Linear size must be positive, got: " + std::to_string(linear_size_));

    if (num_shells_ <= 0)
        throw std::invalid_argument("Number of shells must be positive, got: " + std::to_string(num_shells_));

    if (norm_a_ <= 0 || norm_b_ <= 0)
        throw std::invalid_argument("Lattice constants must be positive, got norm_a: " +
                                    std::to_string(norm_a_) + ", norm_b: " + std::to_string(norm_b_));

    if (translation_.empty() || translation_.size() != 2)
        throw std::runtime_error("Translation vectors not properly initialized");

    const double total_system_x = (std::abs(translation_[0][0]) + std::abs(translation_[1][0])) * linear_size_;
    const double total_system_y = (std::abs(translation_[0][1]) + std::abs(translation_[1][1])) * linear_size_;

    const double max_system_distance = (boundary_conditions_ == BoundaryType::PERIODIC)
                                           ? std::sqrt(total_system_x * total_system_x + total_system_y * total_system_y) / 2.0
                                           : std::sqrt(total_system_x * total_system_x + total_system_y * total_system_y);

    if (static_cast<double>(num_shells_) > max_system_distance)
        throw std::invalid_argument(
            "Search radius (" + std::to_string(num_shells_) +
            ") exceeds system size (" + std::to_string(max_system_distance) +
            "). shells: " + std::to_string(num_shells_) +
            ", min_distance: " + std::to_string(num_shells_) +
            ", linear_size: " + std::to_string(linear_size_));
}

void lattice::Geometry::precomp_indexes()
{
    indexes_.resize(num_atoms_);
    for (int32_t i = 0; i < num_atoms_; ++i)
        indexes_[i] = expand_idx(i);
}

void lattice::Geometry::precomp_coordinates()
{
    coordinates_.resize(num_atoms_);
    for (int32_t i = 0; i < num_atoms_; ++i)
        coordinates_[i] = calculate_coordinate(i);
}

void lattice::Geometry::precomp_neighbors()
{
    neighbors_.resize(num_atoms_);
    for (auto &atom_neighbors : neighbors_)
        atom_neighbors.resize(num_shells_);

    std::vector<std::pair<int32_t, int32_t>> displacements;
    for (int32_t di = -num_shells_; di <= num_shells_; ++di)
        for (int32_t dj = -num_shells_; dj <= num_shells_; ++dj)
            if (di * di + dj * dj <= num_shells_ * num_shells_)
                displacements.emplace_back(di, dj);

    for (int32_t atom_idx = 0; atom_idx < num_atoms_; ++atom_idx)
    {
        std::vector<std::pair<int32_t, double>> atom_neighbors;

        auto [i, j, pos] = indexes_[atom_idx];

        for (const auto &[di, dj] : displacements)
        {
            for (int32_t new_pos = 0; new_pos < num_positions_; ++new_pos)
            {
                int32_t neighbor_idx = collapse_idx(i + di, j + dj, new_pos);

                if (neighbor_idx < 0 || neighbor_idx >= num_atoms_ || neighbor_idx == atom_idx)
                    continue;

                if (boundary_conditions_ == BoundaryType::HARD)
                {
                    auto [ni, nj, npos] = indexes_[neighbor_idx];
                    if (ni < 0 || ni >= linear_size_ || nj < 0 || nj >= linear_size_)
                        continue;
                }

                double distance = calculate_distance(atom_idx, neighbor_idx);
                atom_neighbors.emplace_back(neighbor_idx, distance);
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
            if (std::abs(distance - current_distance) > 1e-5 * current_distance)
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