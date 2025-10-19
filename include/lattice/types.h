#pragma once

namespace lattice
{

    enum class SpinModel
    {
        ISING,
        HEISENBERG
    };

    enum class CrystalType
    {
        RECTANGULAR,
        C_RECTANGULAR,
        TRIANGULAR,
        HONEYCOMB,
        KAGOME,
        LIEB

    };

    enum class BoundaryType
    {
        PERIODIC,
        HARD
    };

} // namespace lattice