#include <filesystem>
#include "../../external/catch2/catch_amalgamated.hpp"
#include "../../include/core/config.h"

TEST_CASE("Test config", "[config]")
{
    Config config;

    auto create_test_config = [](const std::string &path)
    {
        json test_config =
            {
                {"lattice", {{"crystal_type", "rectangular"}, {"system_size", 10}, {"lattice_constant_a", 1.0}, {"lattice_constant_b", 1.0}, {"boundary_type", "periodic"}}},
                {"physical", {{"exchange_constants", {1.0}}, {"external_magnetic_field", {0.0, 0.0, 0.0}}, {"temperature", 10.0}}}};

        std::ofstream file(path);
        file << test_config.dump(4);
        file.close();
    };

    auto create_invalid_config = [](const std::string &path)
    {
        std::ofstream file(path);
        file << "{ invalid json }";
    };

    SECTION("Loading configuration")
    {
        const std::string config_path = "test_config.json";
        create_test_config(config_path);
        config.load(config_path);

        REQUIRE(config.get<int32_t>("lattice.system_size") == 10);
        REQUIRE(config.get<double>("physical.temperature") == 10.0);
        REQUIRE(config.get<std::array<double, 3>>("physical.external_magnetic_field") == std::array<double, 3>{0.0, 0.0, 0.0});

        std::filesystem::remove(config_path);
    }

    SECTION("get/set methods")
    {
        const int32_t TEST_INT_VALUE = 42;
        const double TEST_DOUBLE_VALUE = 3.14;
        const std::string TEST_STRING_VALUE = "Oganechan";
        const double TEST_ZERO_VALUE = 0.0;
        const std::string TEST_EMPTY_STRING = "";

        config.set<int32_t>("test.int_value", TEST_INT_VALUE);
        config.set<double>("test.double_value", TEST_DOUBLE_VALUE);
        config.set<std::string>("test.string_value", TEST_STRING_VALUE);
        config.set<double>("test.zero", 0.0);
        config.set<std::string>("test.empty", "");

        REQUIRE(config.get<int32_t>("test.int_value") == TEST_INT_VALUE);
        REQUIRE(config.get<double>("test.double_value") == TEST_DOUBLE_VALUE);
        REQUIRE(config.get<std::string>("test.string_value") == TEST_STRING_VALUE);
        REQUIRE(config.get<double>("test.TEST_ZERO_VALUE") == TEST_ZERO_VALUE);
        REQUIRE(config.get<std::string>("TEST_EMPTY_STRING") == TEST_EMPTY_STRING);
    }

    SECTION("Non-existent keys, default values")
    {
        const int32_t DEFAULT_INT = 42;
        const double DEFAULT_DOUBLE = 3.14;
        const std::string DEFAULT_STRING = "Sava";

        REQUIRE(config.get<int32_t>("test.non_exist") == 0);
        REQUIRE(config.get<int32_t>("test.non_exist", DEFAULT_INT) == DEFAULT_INT);
        REQUIRE(config.get<double>("test.non_exist_double", DEFAULT_DOUBLE) == DEFAULT_DOUBLE);
        REQUIRE(config.get<std::string>("test.non_exist_string", DEFAULT_STRING) == DEFAULT_STRING);
    }

    SECTION("Error handling")
    {
        const std::string invalid_config_path = "invalid_config.json";
        create_invalid_config(invalid_config_path);

        REQUIRE_THROWS_AS(config.load("outside_of_being.json"), std::runtime_error);
        REQUIRE_THROWS_AS(config.load(invalid_config_path), json::exception);

        std::filesystem::remove(invalid_config_path);
    }
}