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

    inline void set_temperature(double temperature) {
        temperature_ = temperature;
        beta_ = 1.0 / temperature_;
    }

  private:
    lattice::Atoms &atoms_;
    physics::Calculator calculator_;

    double temperature_;
    double beta_;

    inline bool accept_change(double energy_difference) const {
        if (energy_difference <= 0)
            return true;

        const double random_value = Random::uniform_real<double>();
        const double threshold = std::exp(-energy_difference * beta_);

        return random_value <= threshold;
    }
};

} // namespace algorithm

#endif // METROPOLIS_HPP