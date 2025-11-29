#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include "../lattice/atoms.hpp"
#include "../physics/calculator.hpp"

namespace algorithm {

class Metropolis {
  public:
    Metropolis(lattice::Atoms &atoms, const physics::Calculator &calculator);

    void step();
    void sweep();

  private:
    lattice::Atoms &atoms_;
    const physics::Calculator &calculator_;

    bool accept_change(double energy_difference) const;
};

} // namespace algorithm

#endif // METROPOLIS_HPP