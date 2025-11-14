#pragma once

#include <chrono>
#include <omp.h>
#include <random>
#include <thread>

class Random {
  public:
    static std::mt19937 &get_rng() {
        static thread_local std::random_device rd;
        static thread_local std::mt19937 rng(rd());
        return rng;
    }

    static void initialize(uint32_t seed) { get_rng().seed(seed); }

    static void initialize_thread_based() {
        auto seed = std::hash<std::thread::id>{}(std::this_thread::get_id()) ^
                    std::chrono::steady_clock::now().time_since_epoch().count();
        get_rng().seed(seed);
    }

    template <typename T> static T uniform_real(T min = 0.0, T max = 1.0) {
        std::uniform_real_distribution<T> dist(min, max);
        return dist(get_rng());
    }

    template <typename T> static T uniform_int(T min, T max) {
        std::uniform_int_distribution<T> dist(min, max);
        return dist(get_rng());
    }

    static bool bernoulli(double p = 0.5) {
        std::bernoulli_distribution dist(p);
        return dist(get_rng());
    }
};