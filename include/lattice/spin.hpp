#ifndef SPIN_HPP
#define SPIN_HPP

#include "../utils/random.hpp"
#include "types.hpp"
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>

namespace lattice {

class BaseSpin {
  public:
    virtual void randomize() = 0;
    virtual void replace(const BaseSpin &other) = 0;
    virtual void replace(const std::array<double, 3> &other) = 0;
    virtual std::array<double, 3> get_components() const = 0;
    virtual double dot_product(const BaseSpin &other) const = 0;
    virtual double dot_product(const std::array<double, 3> &other) const = 0;

    virtual ~BaseSpin() = default;
};

class IsingSpin : public BaseSpin {
  public:
    void randomize() override { value_ = Random::bernoulli(0.5) ? 1 : -1; }
    void replace(const BaseSpin &other) override {
        value_ = dynamic_cast<const IsingSpin &>(other).value_;
    }
    void replace(const std::array<double, 3> &other) override {
        value_ = other[2] > 0 ? 1 : -1;
    }
    std::array<double, 3> get_components() const override {
        return {0.0, 0.0, double(value_)};
    }
    double dot_product(const BaseSpin &other) const override {
        return value_ * dynamic_cast<const IsingSpin &>(other).value_;
    }
    double dot_product(const std::array<double, 3> &other) const override {
        return value_ * other[2];
    }

  private:
    int32_t value_;
};

class XYSpin : public BaseSpin {
  public:
    void randomize() override { phi_ = Random::uniform_real(0.0, 2.0 * M_PI); }
    void replace(const BaseSpin &other) override {
        phi_ = dynamic_cast<const XYSpin &>(other).phi_;
    }
    void replace(const std::array<double, 3> &other) override {
        phi_ = std::atan2(other[1], other[0]);
    }
    std::array<double, 3> get_components() const override {
        return {std::cos(phi_), std::sin(phi_), 0.0};
    }
    double dot_product(const BaseSpin &other) const override {
        return std::cos(phi_ - dynamic_cast<const XYSpin &>(other).phi_);
    }
    double dot_product(const std::array<double, 3> &other) const override {
        return std::cos(phi_) * other[0] + std::sin(phi_) * other[1];
    }

  private:
    double phi_;
};

class HeisenbergSpin : public BaseSpin {
  public:
    void randomize() override {
        phi_ = Random::uniform_real(0.0, 2.0 * M_PI);
        theta_ = Random::uniform_real(0.0, M_PI);
    }
    void replace(const BaseSpin &other) override {
        const auto &other_heisenberg =
            dynamic_cast<const HeisenbergSpin &>(other);
        phi_ = other_heisenberg.phi_;
        theta_ = other_heisenberg.theta_;
    }
    void replace(const std::array<double, 3> &other) override {
        phi_ = std::atan2(other[1], other[0]);
        theta_ = std::atan2(
            std::sqrt(other[0] * other[0] + other[1] * other[1]), other[2]);
    }
    std::array<double, 3> get_components() const override {
        const double sin_theta = std::sin(theta_);

        return {sin_theta * std::cos(phi_), sin_theta * std::sin(phi_),
                std::cos(theta_)};
    }
    double dot_product(const BaseSpin &other) const override {
        const auto spin1 = get_components();
        const auto spin2 =
            dynamic_cast<const HeisenbergSpin &>(other).get_components();

        return spin1[0] * spin2[0] + spin1[1] * spin2[1] + spin1[2] * spin2[2];
    }
    double dot_product(const std::array<double, 3> &other) const override {
        const auto components = get_components();
        return components[0] * other[0] + components[1] * other[1] +
               components[2] * other[2];
    }

  private:
    double phi_, theta_;
};

class SpinFactory {
  public:
    static std::unique_ptr<BaseSpin> create(SpinModel model) {
        switch (model) {
        case SpinModel::ISING:
            return std::make_unique<IsingSpin>();
        case SpinModel::XY:
            return std::make_unique<XYSpin>();
        case SpinModel::HEISENBERG:
            return std::make_unique<HeisenbergSpin>();
        default:
            throw std::invalid_argument("Unknown spin model");
        }
    }
};

} // namespace lattice

#endif // SPIN_HPP