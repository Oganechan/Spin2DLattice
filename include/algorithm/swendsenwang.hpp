#ifndef SWENDSENWANG_HPP
#define SWENDSENWANG_HPP

#include "../lattice/atoms.hpp"
#include "../physics/calculator.hpp"
#include <cstdint>
#include <vector>

namespace algorithm {

class Swendsenwang {
  public:
    explicit Swendsenwang(lattice::Atoms &atoms,
                          const physics::Calculator &calculator);

    void step();
    void sweep();
    void sweep(int32_t step_count);

  private:
    lattice::Atoms &atoms_;
    const physics::Calculator &calculator_;

    std::vector<int32_t> parent_;
    std::vector<int32_t> rank_;

    int32_t find(int32_t atom_id);
    void unite(int32_t a, int32_t b);
    void build_clusters(const std::array<double, 3> &random_dir);
    void flip_clusters(const std::array<double, 3> &random_dir);
};

} // namespace algorithm

#endif // SWENDSENWANG_HPP