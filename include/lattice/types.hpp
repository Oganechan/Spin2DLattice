#ifndef TYPES_HPP
#define TYPES_HPP

namespace lattice {

enum class CrystalType {
    RECTANGULAR,
    C_RECTANGULAR,
    TRIANGULAR,
    HONEYCOMB,
    KAGOME,
    LIEB

};

enum class BoundaryType { PERIODIC, HARD };

} // namespace lattice

#endif // TYPES_HPP