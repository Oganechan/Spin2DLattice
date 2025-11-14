#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include "../lattice/atoms.hpp"
#include "../physics/calculator.hpp"

namespace algorithm {

class Metropolis {
  public:
    Metropolis(lattice::Atoms &atoms);

    void step();
    void sweep();

  private:
    lattice::Atoms &atoms_;
    physics::Calculator calculator_;

    double temperature_;
    double beta_;

    inline bool accept_change(double energy_difference) const {
        return energy_difference <= 0.0 ||
               Random::uniform_real<double>() <=
                   std::exp(-energy_difference * beta_);
    }
};

} // namespace algorithm

#endif // METROPOLIS_HPP