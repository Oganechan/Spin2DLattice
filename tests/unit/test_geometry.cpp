#include "../../external/catch2/catch_amalgamated.hpp"
#include "../../include/lattice/geometry.hpp"

TEST_CASE("Test geometry", "[geometry]") {
    Config config;

    SECTION("Correct config parsing") {
        config.set<int32_t>("lattice.system_size", 5);
        config.set<double>("lattice.lattice_constant_a", 1.0);
        config.set<double>("lattice.lattice_constant_b", 1.0);
        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_type", "periodic");
        config.set<std::vector<double>>("physical.exchange_constants",
                                        {1.0, 0.25});

        REQUIRE_NOTHROW(lattice::Geometry(config));

        lattice::Geometry geometry(config);
        REQUIRE(geometry.get_system_size() == 5);
        REQUIRE(geometry.get_shell_count() == 2);
        REQUIRE(geometry.get_basis_count() == 1);
        REQUIRE(geometry.get_atom_count() == 25);
        REQUIRE(geometry.get_lattice_constant_a() == 1.0);
        REQUIRE(geometry.get_lattice_constant_b() == 1.0);
        REQUIRE(geometry.get_crystal_type() ==
                lattice::CrystalType::RECTANGULAR);
        REQUIRE(geometry.get_boundary_type() ==
                lattice::BoundaryType::PERIODIC);

        config.set<std::string>("lattice.crystal_type", "invalid_type");
        REQUIRE_THROWS_AS(lattice::Geometry(config), std::invalid_argument);

        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_type", "invalid_boundary");
        REQUIRE_THROWS_AS(lattice::Geometry(config), std::invalid_argument);
    }

    SECTION("Initialize of different crystal lattice") {
        config.set<int32_t>("lattice.system_size", 3);
        config.set<double>("lattice.lattice_constant_a", 1.0);
        config.set<double>("lattice.lattice_constant_b", 1.0);
        config.set<std::string>("lattice.boundary_type", "periodic");
        config.set<std::vector<double>>("physical.exchange_constants", {1.0});

        config.set<std::string>("lattice.crystal_type", "rectangular");
        lattice::Geometry rectangular_geometry(config);
        REQUIRE(rectangular_geometry.get_atom_count() == 3 * 3 * 1);

        config.set<std::string>("lattice.crystal_type", "c_rectangular");
        lattice::Geometry c_rectangular_geometry(config);
        REQUIRE(c_rectangular_geometry.get_atom_count() == 3 * 3 * 2);

        config.set<std::string>("lattice.crystal_type", "triangular");
        lattice::Geometry triangular_geometry(config);
        REQUIRE(triangular_geometry.get_atom_count() == 3 * 3 * 1);

        config.set<std::string>("lattice.crystal_type", "honeycomb");
        lattice::Geometry honeycomb_geometry(config);
        REQUIRE(honeycomb_geometry.get_atom_count() == 3 * 3 * 2);

        config.set<std::string>("lattice.crystal_type", "kagome");
        lattice::Geometry kagome_geometry(config);
        REQUIRE(kagome_geometry.get_atom_count() == 3 * 3 * 3);

        config.set<std::string>("lattice.crystal_type", "lieb");
        lattice::Geometry lieb_geometry(config);
        REQUIRE(lieb_geometry.get_atom_count() == 3 * 3 * 3);
    }

    SECTION("Calculating distances within hard boundaries") {
        config.set<int32_t>("lattice.system_size", 3);
        config.set<double>("lattice.lattice_constant_a", 1.0);
        config.set<double>("lattice.lattice_constant_b", 1.0);
        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_type", "hard");
        config.set<std::vector<double>>("physical.exchange_constants", {1.0});

        lattice::Geometry geometry(config);

        REQUIRE(geometry.get_distance(0, 0) == 0.0);

        const double EXPECTED_DISTANCE = 1.0;
        REQUIRE(geometry.get_distance(4, 1) ==
                Catch::Approx(EXPECTED_DISTANCE).margin(1e-10));
        REQUIRE(geometry.get_distance(4, 3) ==
                Catch::Approx(EXPECTED_DISTANCE).margin(1e-10));
        REQUIRE(geometry.get_distance(4, 5) ==
                Catch::Approx(EXPECTED_DISTANCE).margin(1e-10));
        REQUIRE(geometry.get_distance(4, 7) ==
                Catch::Approx(EXPECTED_DISTANCE).margin(1e-10));

        REQUIRE(geometry.get_distance(0, 1) ==
                Catch::Approx(geometry.get_distance(1, 0)).margin(1e-10));

        const double DIAGONAL_DISTANCE = std::sqrt(2.0);
        REQUIRE(geometry.get_distance(4, 0) ==
                Catch::Approx(DIAGONAL_DISTANCE).margin(1e-10));
        REQUIRE(geometry.get_distance(4, 2) ==
                Catch::Approx(DIAGONAL_DISTANCE).margin(1e-10));
        REQUIRE(geometry.get_distance(4, 6) ==
                Catch::Approx(DIAGONAL_DISTANCE).margin(1e-10));
        REQUIRE(geometry.get_distance(4, 8) ==
                Catch::Approx(DIAGONAL_DISTANCE).margin(1e-10));

        REQUIRE_THROWS_AS(geometry.get_distance(0, -1), std::out_of_range);
        REQUIRE_THROWS_AS(geometry.get_distance(0, 100), std::out_of_range);
    }

    SECTION("Calculating distances within periodic boundaries") {
        config.set<int32_t>("lattice.system_size", 3);
        config.set<double>("lattice.lattice_constant_a", 1.0);
        config.set<double>("lattice.lattice_constant_b", 1.0);
        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_type", "periodic");
        config.set<std::vector<double>>("physical.exchange_constants", {1.0});

        lattice::Geometry geometry(config);

        const double EXPECTED_DISTANCE = 1.0;
        REQUIRE(geometry.get_distance(0, 2) ==
                Catch::Approx(EXPECTED_DISTANCE).margin(1e-10));
        REQUIRE(geometry.get_distance(0, 6) ==
                Catch::Approx(EXPECTED_DISTANCE).margin(1e-10));

        const double DIAGONAL_DISTANCE = std::sqrt(2.0);
        REQUIRE(geometry.get_distance(0, 8) ==
                Catch::Approx(DIAGONAL_DISTANCE).margin(1e-10));
    }

    SECTION("Neighbor table for different types of lattice") {
        config.set<int32_t>("lattice.system_size", 3);
        config.set<double>("lattice.lattice_constant_a", 1.0);
        config.set<double>("lattice.lattice_constant_b", 1.0);
        config.set<std::string>("lattice.boundary_type", "periodic");
        config.set<std::vector<double>>("physical.exchange_constants",
                                        {1.0, 0.25});

        config.set<std::string>("lattice.crystal_type", "rectangular");
        lattice::Geometry rectangular_geometry(config);
        const auto &rectangular_neighbors =
            rectangular_geometry.get_neighbor_table();
        REQUIRE(rectangular_neighbors.size() == 9);
        REQUIRE(rectangular_neighbors[0].size() == 2);
        REQUIRE(rectangular_neighbors[4][1].size() == 4);

        config.set<std::string>("lattice.crystal_type", "honeycomb");
        lattice::Geometry honeycomb_geometry(config);
        const auto &honeycomb_neighbors =
            honeycomb_geometry.get_neighbor_table();
        for (int32_t atom_id = 0; atom_id < honeycomb_geometry.get_atom_count();
             ++atom_id) {
            const auto &atom_neighbors = honeycomb_neighbors[atom_id];
            REQUIRE(atom_neighbors.size() ==
                    honeycomb_geometry.get_shell_count());

            for (int32_t shell = 0;
                 shell < honeycomb_geometry.get_shell_count(); ++shell)
                for (int32_t neighbor_id : atom_neighbors[shell]) {
                    REQUIRE(neighbor_id >= 0);
                    REQUIRE(neighbor_id < honeycomb_geometry.get_atom_count());
                    REQUIRE(neighbor_id != atom_id);
                }
        }

        for (int32_t atom_id = 0; atom_id < honeycomb_geometry.get_atom_count();
             ++atom_id) {
            auto [cell_i, cell_j, atom_in_cell_id] =
                honeycomb_geometry.get_cell_index(atom_id);
            int32_t reconstructed_id =
                honeycomb_geometry.get_atom_id(cell_i, cell_j, atom_in_cell_id);
            REQUIRE(reconstructed_id == atom_id);
        }

        const auto &positions = honeycomb_geometry.get_atom_positions();
        REQUIRE(positions.size() == honeycomb_geometry.get_atom_count());

        for (int32_t atom_id = 0; atom_id < honeycomb_geometry.get_atom_count();
             ++atom_id) {
            auto computed_position =
                honeycomb_geometry.get_atom_position(atom_id);
            REQUIRE(computed_position[0] ==
                    Catch::Approx(positions[atom_id][0]).margin(1e-10));
            REQUIRE(computed_position[1] ==
                    Catch::Approx(positions[atom_id][1]).margin(1e-10));
        }
    }
}