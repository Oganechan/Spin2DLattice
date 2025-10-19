#include "tests.h"

using namespace lattice;

void test_atoms_integration()
{
    print_test_header("Atoms Integration Test");

    Config config;
    fill_test_config(config);

    Atoms atoms(config);

    std::cout << "Model type: " << static_cast<int32_t>(atoms.get_type()) << std::endl;
    std::cout << "Geometry atoms: " << atoms.geometry().get_num_atoms() << std::endl;

    atoms.random_initialize();
    std::cout << "Initialization completed" << std::endl;

    atoms.set_magnetic(1, false);
    std::cout << "Atom 1 magnetic: " << (atoms.is_magnetic(1) ? "Yes" : "No") << std::endl;

    print_success("Atoms test PASSED");
}

int main()
{
    std::cout << "Starting Atoms Test..." << std::endl
              << std::endl;

    try
    {
        test_atoms_integration();
    }
    catch (const std::exception &e)
    {
        std::cerr << "TEST FAILED: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}