#include "../../external/catch2/catch_amalgamated.hpp"
#include "../../include/lattice/atoms.h"

TEST_CASE("Test atoms", "[atoms]")
{
    Config config;
    config.set<int32_t>("lattice.system_size", 4);
    config.set<double>("lattice.lattice_constant_a", 1.0);
    config.set<double>("lattice.lattice_constant_b", 1.0);
    config.set<std::string>("lattice.crystal_type", "rectangular");
    config.set<std::string>("lattice.boundary_type", "periodic");
    config.set<std::vector<double>>("physical.exchange_constants", {1.0, 0.5});

    SECTION("Correct initialize")
    {
        lattice::Atoms atoms(config);

        REQUIRE(atoms.get_geometry().get_atom_count() == 16);
        REQUIRE(atoms.get_magnetic_count() == 16);
        REQUIRE(atoms.get_defect_count() == 0);
    }

    SECTION("Magnetic mask tests")
    {
        lattice::Atoms atoms(config);

        atoms.set_magnetic_state(0, false);
        atoms.set_magnetic_state(1, true);

        REQUIRE_FALSE(atoms.get_magnetic_state(0));
        REQUIRE(atoms.get_magnetic_state(1));

        REQUIRE_THROWS_AS(atoms.set_magnetic_state(-1, true), std::out_of_range);
        REQUIRE_THROWS_AS(atoms.set_magnetic_state(100, true), std::out_of_range);
        REQUIRE_THROWS_AS(atoms.get_magnetic_state(-1), std::out_of_range);
        REQUIRE_THROWS_AS(atoms.get_magnetic_state(100), std::out_of_range);

        auto mask = atoms.get_magnetic_mask();
        REQUIRE(mask.size() == 16);
        REQUIRE_FALSE(mask[0]);
        REQUIRE(mask[1]);
    }

    SECTION("Random defect generation")
    {
        lattice::Atoms atoms(config);

        REQUIRE_NOTHROW(atoms.set_random_defects(0.25));
        REQUIRE_NOTHROW(atoms.set_random_defects(0.0));
        REQUIRE_NOTHROW(atoms.set_random_defects(1.0));

        REQUIRE_THROWS_AS(atoms.set_random_defects(-0.1), std::invalid_argument);
        REQUIRE_THROWS_AS(atoms.set_random_defects(1.1), std::invalid_argument);

        atoms.set_random_defects(0.25);
        const int32_t EXPECTED_DEFECTS = static_cast<int32_t>(0.25 * 16);
        REQUIRE(atoms.get_defect_count() == EXPECTED_DEFECTS);
        REQUIRE(atoms.get_magnetic_count() == 16 - EXPECTED_DEFECTS);
        REQUIRE(atoms.get_magnetic_count() + atoms.get_defect_count() == 16);
    }

    SECTION("Spin management")
    {
        lattice::Atoms atoms(config);

        std::array<double, 3> spin = {0.6, 0.8, 0.0};
        atoms.set_spin(0, spin);

        auto retrieved_spin = atoms.get_spin(0);
        double norm = std::sqrt(retrieved_spin[0] * retrieved_spin[0] +
                                retrieved_spin[1] * retrieved_spin[1] +
                                retrieved_spin[2] * retrieved_spin[2]);
        REQUIRE(norm == Catch::Approx(1.0).margin(1e-10));

        atoms.set_magnetic_state(1, false);
        REQUIRE_THROWS_AS(atoms.set_spin(1, spin), std::invalid_argument);
        REQUIRE_THROWS_AS(atoms.get_spin(1), std::invalid_argument);

        REQUIRE_THROWS_AS(atoms.set_spin(-1, spin), std::out_of_range);
        REQUIRE_THROWS_AS(atoms.get_spin(-1), std::out_of_range);

        std::array<double, 3> unnormalized_spin = {2.0, 3.0, 4.0};
        atoms.set_spin(2, unnormalized_spin);
        auto normalized_spin = atoms.get_spin(2);
        double normalized_norm = std::sqrt(normalized_spin[0] * normalized_spin[0] +
                                           normalized_spin[1] * normalized_spin[1] +
                                           normalized_spin[2] * normalized_spin[2]);
        REQUIRE(normalized_norm == Catch::Approx(1.0).margin(1e-10));
    }

    SECTION("Spin initialization methods")
    {
        lattice::Atoms atoms(config);

        atoms.initialize_random();
        bool found_different_spins = false;
        auto first_spin = atoms.get_spin(0);
        for (int32_t i = 1; i < 16; ++i)
        {
            auto current_spin = atoms.get_spin(i);
            if (std::abs(current_spin[0] - first_spin[0]) > 1e-10 ||
                std::abs(current_spin[1] - first_spin[1]) > 1e-10 ||
                std::abs(current_spin[2] - first_spin[2]) > 1e-10)
            {
                found_different_spins = true;
                break;
            }
        }
        REQUIRE(found_different_spins);

        atoms.initialize_ferromagnetic();
        for (int32_t i = 0; i < 16; ++i)
        {
            auto spin = atoms.get_spin(i);
            REQUIRE(spin[0] == Catch::Approx(0.0).margin(1e-10));
            REQUIRE(spin[1] == Catch::Approx(0.0).margin(1e-10));
            REQUIRE(spin[2] == Catch::Approx(1.0).margin(1e-10));
        }

        atoms.initialize_antiferromagnetic();
        int32_t up_count = 0, down_count = 0;
        for (int32_t i = 0; i < 16; ++i)
        {
            auto spin = atoms.get_spin(i);
            REQUIRE(spin[0] == Catch::Approx(0.0).margin(1e-10));
            REQUIRE(spin[1] == Catch::Approx(0.0).margin(1e-10));
            if (spin[2] > 0.5)
                up_count++;
            else if (spin[2] < -0.5)
                down_count++;
        }
        REQUIRE(up_count > 0);
        REQUIRE(down_count > 0);
        REQUIRE(up_count + down_count == 16);
    }
}