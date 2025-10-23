#pragma once

#include "../lattice/atoms.h"
#include "../physics/calculator.h"

namespace algorithm
{

    class Metropolis
    {
    public:
        Metropolis(lattice::Atoms &atoms, double temperature);

        void step();
        void sweep();

    private:
        lattice::Atoms &atoms_;
        physics::Calculator calculator_;

        double temperature_;
        double beta_;

        int32_t select_random_magnetic_atom() const;
        inline bool accept_change(double energy_difference) const
        {
            return energy_difference <= 0.0 || Random::uniform_real<double>() <= std::exp(-energy_difference * beta_);
        }
    };

} // namespace algorithm