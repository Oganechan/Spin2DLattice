#ifndef SWENDSENWANG_HPP
#define SWENDSENWANG_HPP

#include "../lattice/atoms.hpp"
#include "../physics/calculator.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace algorithm {

class ISwendsenWang {
  public:
    virtual ~ISwendsenWang() = default;
    virtual void step() = 0;
    virtual void sweep() = 0;
    virtual void sweep(int32_t step_count) = 0;
};

class SwendsenWangBase : public ISwendsenWang {
  protected:
    lattice::Atoms &atoms_;
    const physics::Calculator &calculator_;

  public:
    SwendsenWangBase(lattice::Atoms &atoms,
                     const physics::Calculator &calculator)
        : atoms_(atoms), calculator_(calculator) {}

    virtual ~SwendsenWangBase() = default;

    void sweep() override {
        int32_t step_count = atoms_.get_magnetic_count() / 10;
        while (step_count-- > 0)
            step();
    }

    void sweep(int32_t step_count) override {
        while (step_count-- > 0)
            step();
    }
};

class SwendsenWangIsing : public SwendsenWangBase {
  public:
    SwendsenWangIsing(lattice::Atoms &atoms,
                      const physics::Calculator calculator)
        : SwendsenWangBase(atoms, calculator) {}

    void step() override {
        const auto dir = generate_random_direction();
        build_clusters(dir);
        flip_clusters(dir);
    }

  private:
    std::vector<int32_t> parent_;
    std::vector<int32_t> rank_;

    std::array<double, 3> generate_random_direction() {
        return {0.0, 0.0, Random::bernoulli(0.5) ? 1.0 : -1.0};
    }

    int32_t find(int32_t atom_id) {
        while (parent_[atom_id] != atom_id) {
            parent_[atom_id] = parent_[parent_[atom_id]];
            atom_id = parent_[atom_id];
        }

        return atom_id;
    }

    void unite(int32_t a, int32_t b) {
        int32_t first_parent = find(a);
        int32_t second_parent = find(b);

        if (first_parent == second_parent)
            return;

        if (rank_[first_parent] < rank_[second_parent])
            parent_[first_parent] = second_parent;
        else {
            parent_[second_parent] = first_parent;
            if (rank_[first_parent] == rank_[second_parent])
                rank_[first_parent]++;
        }
    }

    void build_clusters(const std::array<double, 3> &random_dir) {
        const int32_t atom_count = atoms_.get_geometry().get_atom_count();

        parent_.resize(atom_count);
        rank_.resize(atom_count, 0);
        for (int32_t i = 0; i < atom_count; ++i)
            parent_[i] = i;

        std::vector<double> projections(atom_count, 0.0);
        for (int32_t atom_id = 0; atom_id < atom_count; ++atom_id)
            if (atoms_.get_magnetic_state(atom_id)) {
                const auto &spin = atoms_.get_spin(atom_id).get_components();
                projections[atom_id] = spin[0] * random_dir[0] +
                                       spin[1] * random_dir[1] +
                                       spin[2] * random_dir[2];
            }

        for (int32_t atom_id = 0; atom_id < atom_count; ++atom_id) {
            for (int32_t shell = 0;
                 shell < atoms_.get_geometry().get_shell_count(); ++shell) {
                const double J = calculator_.get_exchange_constant(shell);

                for (int32_t neighbor_id :
                     atoms_.get_geometry()
                         .get_neighbor_table()[atom_id][shell]) {
                    if (neighbor_id <= atom_id)
                        continue;
                    if (!atoms_.get_magnetic_state(neighbor_id))
                        continue;

                    double p_bond =
                        1.0 - std::exp(-2.0 * J * projections[atom_id] *
                                       projections[neighbor_id] /
                                       calculator_.get_temperature());
                    p_bond = std::max(0.0, std::min(1.0, p_bond));

                    if (Random::bernoulli(p_bond))
                        unite(atom_id, neighbor_id);
                }
            }
        }
    }

    void flip_clusters(const std::array<double, 3> &random_dir) {
        std::unordered_map<int32_t, std::vector<int32_t>> cluster_map;

        for (int32_t atom_id : atoms_.get_magnetic_atoms()) {
            int32_t root = find(atom_id);
            cluster_map[root].push_back(atom_id);
        }

        for (auto &[root, cluster_atoms] : cluster_map) {
            if (Random::bernoulli()) {
                for (int32_t atom_id : cluster_atoms)
                    flip_spin(atom_id, random_dir);
            }
        }
    }

    void flip_spin(int32_t atom_id, const std::array<double, 3> &dir) {
        const auto &spin = atoms_.get_spin(atom_id);
        const auto components = spin.get_components();
        std::array<double, 3> new_components = {0.0, 0.0, -components[2]};

        std::unique_ptr<lattice::BaseSpin> new_spin =
            atoms_.generate_random_spin();
        new_spin->replace(new_components);
        atoms_.set_spin(atom_id, *new_spin);
    }
};

class SwendsenWangXY : public SwendsenWangBase {
  public:
    SwendsenWangXY(lattice::Atoms &atoms, const physics::Calculator &calculator)
        : SwendsenWangBase(atoms, calculator) {}

    void step() override {
        const auto dir = generate_random_direction();
        build_clusters(dir);
        flip_clusters(dir);
    }

  private:
    std::vector<int32_t> parent_;
    std::vector<int32_t> rank_;

    std::array<double, 3> generate_random_direction() {
        double phi = Random::uniform_real(0.0, 2.0 * M_PI);
        return {std::cos(phi), std::sin(phi), 0.0};
    }

    int32_t find(int32_t atom_id) {
        while (parent_[atom_id] != atom_id) {
            parent_[atom_id] = parent_[parent_[atom_id]];
            atom_id = parent_[atom_id];
        }

        return atom_id;
    }

    void unite(int32_t a, int32_t b) {
        int32_t first_parent = find(a);
        int32_t second_parent = find(b);

        if (first_parent == second_parent)
            return;

        if (rank_[first_parent] < rank_[second_parent])
            parent_[first_parent] = second_parent;
        else {
            parent_[second_parent] = first_parent;
            if (rank_[first_parent] == rank_[second_parent])
                rank_[first_parent]++;
        }
    }

    void build_clusters(const std::array<double, 3> &random_dir) {
        const int32_t atom_count = atoms_.get_geometry().get_atom_count();

        parent_.resize(atom_count);
        rank_.resize(atom_count, 0);
        for (int32_t i = 0; i < atom_count; ++i)
            parent_[i] = i;

        std::vector<double> projections(atom_count, 0.0);
        for (int32_t atom_id = 0; atom_id < atom_count; ++atom_id)
            if (atoms_.get_magnetic_state(atom_id)) {
                const auto &spin = atoms_.get_spin(atom_id).get_components();
                projections[atom_id] = spin[0] * random_dir[0] +
                                       spin[1] * random_dir[1] +
                                       spin[2] * random_dir[2];
            }

        for (int32_t atom_id = 0; atom_id < atom_count; ++atom_id) {
            for (int32_t shell = 0;
                 shell < atoms_.get_geometry().get_shell_count(); ++shell) {
                const double J = calculator_.get_exchange_constant(shell);

                for (int32_t neighbor_id :
                     atoms_.get_geometry()
                         .get_neighbor_table()[atom_id][shell]) {
                    if (neighbor_id <= atom_id)
                        continue;
                    if (!atoms_.get_magnetic_state(neighbor_id))
                        continue;

                    double p_bond =
                        1.0 - std::exp(-2.0 * J * projections[atom_id] *
                                       projections[neighbor_id] /
                                       calculator_.get_temperature());
                    p_bond = std::max(0.0, std::min(1.0, p_bond));

                    if (Random::bernoulli(p_bond))
                        unite(atom_id, neighbor_id);
                }
            }
        }
    }

    void flip_clusters(const std::array<double, 3> &random_dir) {
        std::unordered_map<int32_t, std::vector<int32_t>> cluster_map;

        for (int32_t atom_id : atoms_.get_magnetic_atoms()) {
            int32_t root = find(atom_id);
            cluster_map[root].push_back(atom_id);
        }

        for (auto &[root, cluster_atoms] : cluster_map) {
            if (Random::bernoulli()) {
                for (int32_t atom_id : cluster_atoms)
                    flip_spin(atom_id, random_dir);
            }
        }
    }

    void flip_spin(int32_t atom_id, const std::array<double, 3> &dir) {
        const auto &spin = atoms_.get_spin(atom_id);
        const auto components = spin.get_components();

        double dot = components[0] * dir[0] + components[1] * dir[1];
        std::array<double, 3> new_components = {
            components[0] - 2.0 * dot * dir[0],
            components[1] - 2.0 * dot * dir[1], 0.0};

        std::unique_ptr<lattice::BaseSpin> new_spin =
            atoms_.generate_random_spin();
        new_spin->replace(new_components);
        atoms_.set_spin(atom_id, *new_spin);
    }
};

class SwendsenWangHeisenberg : public SwendsenWangBase {
  public:
    SwendsenWangHeisenberg(lattice::Atoms &atoms,
                           const physics::Calculator &calculator)
        : SwendsenWangBase(atoms, calculator) {}

    void step() override {
        const auto dir = generate_random_direction();
        build_clusters(dir);
        flip_clusters(dir);
    }

  private:
    std::vector<int32_t> parent_;
    std::vector<int32_t> rank_;

    std::array<double, 3> generate_random_direction() {
        double phi = Random::uniform_real(0.0, 2.0 * M_PI);
        double theta = Random::uniform_real(0.0, M_PI);
        double sin_theta = std::sin(theta);

        return {std::cos(phi) * sin_theta, std::sin(phi) * sin_theta,
                std::cos(theta)};
    }

    int32_t find(int32_t atom_id) {
        while (parent_[atom_id] != atom_id) {
            parent_[atom_id] = parent_[parent_[atom_id]];
            atom_id = parent_[atom_id];
        }

        return atom_id;
    }

    void unite(int32_t a, int32_t b) {
        int32_t first_parent = find(a);
        int32_t second_parent = find(b);

        if (first_parent == second_parent)
            return;

        if (rank_[first_parent] < rank_[second_parent])
            parent_[first_parent] = second_parent;
        else {
            parent_[second_parent] = first_parent;
            if (rank_[first_parent] == rank_[second_parent])
                rank_[first_parent]++;
        }
    }

    void build_clusters(const std::array<double, 3> &random_dir) {
        const int32_t atom_count = atoms_.get_geometry().get_atom_count();

        parent_.resize(atom_count);
        rank_.resize(atom_count, 0);
        for (int32_t i = 0; i < atom_count; ++i)
            parent_[i] = i;

        std::vector<double> projections(atom_count, 0.0);
        for (int32_t atom_id = 0; atom_id < atom_count; ++atom_id)
            if (atoms_.get_magnetic_state(atom_id)) {
                const auto &spin = atoms_.get_spin(atom_id).get_components();
                projections[atom_id] = spin[0] * random_dir[0] +
                                       spin[1] * random_dir[1] +
                                       spin[2] * random_dir[2];
            }

        for (int32_t atom_id = 0; atom_id < atom_count; ++atom_id) {
            for (int32_t shell = 0;
                 shell < atoms_.get_geometry().get_shell_count(); ++shell) {
                const double J = calculator_.get_exchange_constant(shell);

                for (int32_t neighbor_id :
                     atoms_.get_geometry()
                         .get_neighbor_table()[atom_id][shell]) {
                    if (neighbor_id <= atom_id)
                        continue;
                    if (!atoms_.get_magnetic_state(neighbor_id))
                        continue;

                    double p_bond =
                        1.0 - std::exp(-2.0 * J * projections[atom_id] *
                                       projections[neighbor_id] /
                                       calculator_.get_temperature());
                    p_bond = std::max(0.0, std::min(1.0, p_bond));

                    if (Random::bernoulli(p_bond))
                        unite(atom_id, neighbor_id);
                }
            }
        }
    }

    void flip_clusters(const std::array<double, 3> &random_dir) {
        std::unordered_map<int32_t, std::vector<int32_t>> cluster_map;

        for (int32_t atom_id : atoms_.get_magnetic_atoms()) {
            int32_t root = find(atom_id);
            cluster_map[root].push_back(atom_id);
        }

        for (auto &[root, cluster_atoms] : cluster_map) {
            if (Random::bernoulli()) {
                for (int32_t atom_id : cluster_atoms)
                    flip_spin(atom_id, random_dir);
            }
        }
    }

    void flip_spin(int32_t atom_id, const std::array<double, 3> &dir) {
        const auto &spin = atoms_.get_spin(atom_id);
        const auto components = spin.get_components();

        double dot = components[0] * dir[0] + components[1] * dir[1] +
                     components[2] * dir[2];
        std::array<double, 3> new_components = {
            components[0] - 2.0 * dot * dir[0],
            components[1] - 2.0 * dot * dir[1],
            components[2] - 2.0 * dot * dir[2]};

        std::unique_ptr<lattice::BaseSpin> new_spin =
            atoms_.generate_random_spin();
        new_spin->replace(new_components);
        atoms_.set_spin(atom_id, *new_spin);
    }
};

class SwendsenWangFactory {
  public:
    static std::unique_ptr<ISwendsenWang>
    create(uint32_t dir_count, lattice::Atoms &atoms,
           const physics::Calculator &calculator) {
        switch (dir_count) {
        case 1:
            return std::make_unique<SwendsenWangIsing>(atoms, calculator);
        case 2:
            return std::make_unique<SwendsenWangXY>(atoms, calculator);
        case 3:
            return std::make_unique<SwendsenWangHeisenberg>(atoms, calculator);
        default:
            throw std::runtime_error("Unsupported direction count");
        }
    }
};

} // namespace algorithm

#endif // SWENDSENWANG_HPP