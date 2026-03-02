#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cstdint>
#include <random>
#include <stdexcept>
#include <type_traits>
class Random {
  public:
    Random(const Random &) = delete;
    Random &operator=(const Random &) = delete;
    Random(Random &&) = default;
    Random &operator=(Random &&) = default;

    explicit Random() { reseed(); }
    explicit Random(uint32_t sd) : engine_(sd) {}

    template <typename T> T uniform_int(T min, T max) {
        static_assert(std::is_integral<T>::value);
        std::uniform_int_distribution<T> dist(min, max);
        return dist(engine_);
    }

    template <typename T> T uniform_real(T min, T max) {
        static_assert(std::is_floating_point<T>::value);
        std::uniform_real_distribution<T> dist(min, max);
        return dist(engine_);
    }

    bool bernoulli(double p = 0.5) {
        if (p < 0.0 || p > 1.0)
            throw std::out_of_range("");
        std::bernoulli_distribution dist(p);
        return dist(engine_);
    }

  private:
    std::mt19937 engine_;

    void reseed() {
        std::random_device rd;
        engine_.seed(rd());
    }
    void reseed(uint32_t sd) { engine_.seed(sd); }
};

inline Random &thread_local_random() {
    thread_local Random rng;
    return rng;
}

inline Random &thread_local_random(uint32_t sd) {
    thread_local Random rng(sd);
    return rng;
}

#endif // RANDOM_HPP