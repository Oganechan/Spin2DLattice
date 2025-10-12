#pragma once

#include <random>

class Random
{
public:
    static std::mt19937 &get_rng()
    {
        static thread_local std::random_device rd;
        static thread_local std::mt19937 rng(rd());
        return rng;
    }
};