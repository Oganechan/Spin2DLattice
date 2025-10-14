#pragma once

namespace lattice
{

    enum class SpinModel
    {
        ISING,
        HEISENBERG,
        XY
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

}