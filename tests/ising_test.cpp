#include "tests.h"
#include "../include/lattice/models/ising_model.h"

using namespace lattice;

void test_basis_initialization()
{
    print_test_header("Basis Initialization");

    Config config;
    fill_test_config(config, 3, "ising");

    Atoms atoms(config);
    Geometry geometry = atoms.geometry();

    std::cout << "Num atoms: " << geometry.get_num_atoms() << std::endl
              << "Model Type: " << static_cast<int32_t>(atoms.get_type()) << std::endl;

    const auto &mask = atoms.get_magnetic_mask();
    bool all_magnetic = true;
    for (bool magnetic : mask)
        if (!magnetic)
        {
            all_magnetic = false;
            break;
        }

    std::cout << "All atoms initially magnetic: " << (all_magnetic ? "YES" : "NO") << std::endl;

    print_success("Basis initialization completed");
}

void test_spin_operations()
{
    print_test_header("Ising Spin Operations");

    Config config;
    fill_test_config(config, 2, "ising");

    Atoms atoms(config);
    Geometry geometry = atoms.geometry();

    atoms.set_spin(0, 1);
    atoms.set_spin(1, -1);
    atoms.set_spin(2, 1);
    atoms.set_spin(3, -1);

    std::cout << "Spin values: ";
    for (int32_t i = 0; i < 4; ++i)
    {
        int32_t spin = std::get<int32_t>(atoms.get_spin(i));
        std::cout << "s" << i << "=" << spin << " ";
    }
    std::cout << std::endl;

    print_success("Spin operations test completed");
}

void test_magnetic_defects()
{
    print_test_header("Ising Magnetic Defects");

    Config config;
    fill_test_config(config, 4, "ising");

    Atoms atoms(config);
    Geometry geometry = atoms.geometry();

    atoms.set_magnetic(2, false);
    atoms.set_magnetic(5, false);
    atoms.set_magnetic(8, false);
    atoms.set_magnetic(12, false);

    std::cout << "Magnetic atoms count: " << atoms.count_magnetic_atoms() << std::endl
              << "Defects count: " << atoms.count_defects() << std::endl;

    auto magnetic_indices = atoms.get_magnetic_atom_indices();
    auto defect_indices = atoms.get_defect_indices();

    print_vector(magnetic_indices, "Magnetic indices");
    print_vector(defect_indices, "Defect indices");

    std::cout << "Checking defect behavior:" << std::endl;
    for (int32_t defect_idx : defect_indices)
    {
        bool is_mag = atoms.is_magnetic(defect_idx);
        std::cout << "Atom " << defect_idx << " is magnetic: " << (is_mag ? "YES" : "NO") << std::endl;
    }

    print_success("Magnetic defects test completed");
}

void test_random_defects()
{
    print_test_header("Ising Random Defects");

    Config config;
    fill_test_config(config, 5, "ising");

    Atoms atoms(config);
    Geometry geometry = atoms.geometry();

    std::vector<double> concentrations = {0.0, 0.2, 0.5, 1.0};

    for (double concentration : concentrations)
    {
        atoms.set_random_defects(concentration);

        int32_t magnetic_count = atoms.count_magnetic_atoms();
        int32_t defect_count = atoms.count_defects();
        double actual_concentration = static_cast<double>(defect_count) / geometry.get_num_atoms();

        std::cout << "Requested: " << concentration << " | actual: " << actual_concentration;
        std::cout << " | Magnetic: " << magnetic_count << " | Defects: " << defect_count;
        std::cout << " | Total: " << geometry.get_num_atoms() << std::endl;
    }

    print_success("Random defects test completed");
}

void test_random_initialization()
{
    print_test_header("Ising Random Initialization");

    Config config;
    fill_test_config(config, 4, "ising");

    Atoms atoms(config);
    Geometry geometry = atoms.geometry();

    atoms.random_initialize();
    atoms.set_random_defects(0.2);

    auto magnetic_indices = atoms.get_magnetic_atom_indices();
    int32_t positive_count = 0;
    int32_t negative_count = 0;

    for (int32_t idx : magnetic_indices)
    {
        int32_t spin = std::get<int32_t>(atoms.get_spin(idx));
        if (spin == 1)
            positive_count++;
        else if (spin == -1)
            negative_count++;
    }

    std::cout << "Magnetic atoms: " << magnetic_indices.size() << std::endl
              << "Positive spins: " << positive_count << std::endl
              << "Negative spins: " << negative_count << std::endl
              << "Ratio: " << static_cast<double>(positive_count) / magnetic_indices.size() << std::endl;

    print_success("Random initialization test completed");
}

int main()
{
    std::cout << "STARTING ISING MODEL TEST" << std::endl
              << "=========================" << std::endl;

    try
    {
        test_basis_initialization();
        test_spin_operations();
        test_magnetic_defects();
        test_random_defects();
        test_random_initialization();

        std::cout << std::endl
                  << std::string(50, '=') << std::endl;
        std::cout << "ALL ISING MODEL TESTS COMPLETED SUCCESSFULLY!" << std::endl;
        std::cout << std::string(50, '=') << std::endl;
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "TEST FAILED: " << e.what() << std::endl;
        return 1;
    }
}