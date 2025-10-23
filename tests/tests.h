#pragma once

#include <iostream>
#include "../include/core/simulation.h"

void print_test_header(const std::string &test_name)
{
    std::cout << std::endl
              << std::string(50, '=') << std::endl;

    std::cout << "TEST: " << test_name;

    std::cout << std::endl
              << std::string(50, '=') << std::endl;
}

void print_success(const std::string &message)
{
    std::cout << "✓ " << message << std::endl;
}

void print_failure(const std::string &message)
{
    std::cout << "✗ " << message << std::endl;
}

void print_vector(const std::vector<int32_t> &vec, const std::string &name)
{
    int32_t max_num_args = 10;
    std::cout << name << " [" << vec.size() << "]: ";
    for (int32_t i = 0; i < vec.size() && i < max_num_args; ++i)
    {
        std::cout << vec[i];
        if (i < vec.size() - 1 && i < max_num_args - 1)
            std::cout << ", ";
    }
    if (vec.size() > max_num_args)
        std::cout << "...";
    std::cout << std::endl;
}

void print_vector(const std::vector<bool> &vec, const std::string &name)
{
    int32_t max_num_args = 10;
    std::cout << name << " [" << vec.size() << "]: ";
    for (int32_t i = 0; i < vec.size() && i < max_num_args; ++i)
    {
        std::cout << (vec[i] ? "T" : "F");
        if (i < vec.size() - 1 && i < max_num_args - 1)
            std::cout << ", ";
    }
    if (vec.size() >> max_num_args)
        std::cout << "...";
    std::cout << std::endl;
}

void fill_test_config(Config &config, int32_t system_size = 5)
{

    config.set("lattice.system_size", system_size);
    config.set("lattice.lattice_constant_a", 1.0);
    config.set("lattice.lattice_constant_b", 1.0);
    config.set("lattice.crystal_type", "rectangular");
    config.set("lattice.boundary_type", "periodic");
    config.set("physical.exchange_constants", std::vector<double>{1.0});
}