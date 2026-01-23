#include "swendsenwang.hpp"
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

algorithm::Swendsenwang::Swendsenwang(lattice::Atoms &atoms,
                                      const physics::Calculator &calculator)
    : atoms_(atoms), calculator_(calculator) {}

void algorithm::Swendsenwang::step() {
    const std::array<double, 3> random_dir = atoms_.generate_random_spin();

    build_clusters(random_dir);
    flip_clusters(random_dir);
}

void algorithm::Swendsenwang::sweep() {
    int32_t step_count = atoms_.get_magnetic_count() / 10;

    while (step_count-- > 0)
        step();
}

void algorithm::Swendsenwang::sweep(int32_t step_count) {
    while (step_count-- > 0)
        step();
}

int32_t algorithm::Swendsenwang::find(int32_t atom_id) {
    while (parent_[atom_id] != atom_id) {
        parent_[atom_id] = parent_[parent_[atom_id]];
        atom_id = parent_[atom_id];
    }

    return atom_id;
}

void algorithm::Swendsenwang::unite(int32_t a, int32_t b) {
    int32_t first_parent = find(a);
    int32_t second_parent = find(b);

    if (first_parent == second_parent)
        return;

    if (rank_[a] < rank_[b])
        parent_[first_parent] = second_parent;
    else {
        parent_[second_parent] = first_parent;
        if (rank_[first_parent] == rank_[second_parent])
            rank_[first_parent]++;
    }
}

void algorithm::Swendsenwang::build_clusters(
    const std::array<double, 3> &random_dir) {
    const int32_t atom_count = atoms_.get_geometry().get_atom_count();

    parent_.resize(atom_count);
    rank_.resize(atom_count, 0);
    for (int32_t i = 0; i < atom_count; ++i)
        parent_[i] = i;

    std::vector<double> projections(atom_count);
    for (int32_t atom_id = 0; atom_id < atom_count; ++atom_id) {
        if (atoms_.get_magnetic_state(atom_id)) {
            const auto spin = atoms_.get_spin(atom_id);
            projections[atom_id] = spin[0] * random_dir[0] +
                                   spin[1] * random_dir[1] +
                                   spin[2] * random_dir[2];
        } else
            projections[atom_id] = 0.0;
    }

    for (int32_t atom_id : atoms_.get_magnetic_atoms())
        for (int32_t shell = 0; shell < atoms_.get_geometry().get_shell_count();
             ++shell)
            for (int32_t neighbor_id :
                 atoms_.get_geometry().get_neighbor_table()[atom_id][shell]) {
                if (neighbor_id <= atom_id)
                    continue;
                if (!atoms_.get_magnetic_state(neighbor_id))
                    continue;

                double p_bond =
                    1.0 -
                    std::exp(-2.0 * calculator_.get_exchange_constant(shell) *
                             projections[atom_id] * projections[neighbor_id] /
                             calculator_.get_temperature());
                p_bond = std::max(0.0, std::min(1.0, p_bond));

                if (Random::bernoulli(p_bond))
                    unite(atom_id, neighbor_id);
            }
}

void algorithm::Swendsenwang::flip_clusters(
    const std::array<double, 3> &random_dir) {
    std::unordered_map<int32_t, std::vector<int32_t>> cluster_map;

    for (int32_t atom_id : atoms_.get_magnetic_atoms()) {
        int32_t root = find(atom_id);
        cluster_map[root].push_back(atom_id);
    }

    for (auto &[root, cluster_atoms] : cluster_map) {
        if (Random::bernoulli())
            for (int32_t atom_id : cluster_atoms) {
                const auto &old_spin = atoms_.get_spin(atom_id);

                double projections = old_spin[0] * random_dir[0] +
                                     old_spin[1] * random_dir[1] +
                                     old_spin[2] * random_dir[2];

                std::array<double, 3> new_spin = {
                    old_spin[0] - 2.0 * projections * random_dir[0],
                    old_spin[1] - 2.0 * projections * random_dir[1],
                    old_spin[2] - 2.0 * projections * random_dir[2]};

                atoms_.set_spin(atom_id, new_spin);
            }
    }
}