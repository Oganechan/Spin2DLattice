#ifndef SPIN_HPP
#define SPIN_HPP

#include "../utils/random.hpp"
#include "types.hpp"
#include <array>
#include <cmath>
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

  protected:
    Random &rng() { return thread_local_random(); }
};

class IsingSpin : public BaseSpin {
  public:
    explicit IsingSpin(double spin_value) : spin_value_(spin_value) {}

    void randomize() override {
        value_ = rng().bernoulli(0.5) ? spin_value_ : -spin_value_;
    }
    void replace(const BaseSpin &other) override {
        value_ = dynamic_cast<const IsingSpin &>(other).value_;
    }
    void replace(const std::array<double, 3> &other) override {
        value_ =
            other[2] > 0 ? spin_value_ : (other[2] < 0 ? -spin_value_ : value_);
    }
    std::array<double, 3> get_components() const override {
        return {0.0, 0.0, value_};
    }
    double dot_product(const BaseSpin &other) const override {
        return value_ * dynamic_cast<const IsingSpin &>(other).value_;
    }
    double dot_product(const std::array<double, 3> &other) const override {
        return value_ * other[2];
    }

  private:
    double value_;
    const double spin_value_;
};

class XYSpin : public BaseSpin {
  public:
    explicit XYSpin(double spin_value) : spin_value_(spin_value) {}

    void randomize() override { phi_ = rng().uniform_real(0.0, 2.0 * M_PI); }
    void replace(const BaseSpin &other) override {
        phi_ = dynamic_cast<const XYSpin &>(other).phi_;
    }
    void replace(const std::array<double, 3> &other) override {
        phi_ = std::atan2(other[1], other[0]);
    }
    std::array<double, 3> get_components() const override {
        return {spin_value_ * std::cos(phi_), spin_value_ * std::sin(phi_),
                0.0};
    }
    double dot_product(const BaseSpin &other) const override {
        return spin_value_ *
               std::cos(phi_ - dynamic_cast<const XYSpin &>(other).phi_);
    }
    double dot_product(const std::array<double, 3> &other) const override {
        return spin_value_ * std::cos(phi_) * other[0] +
               spin_value_ * std::sin(phi_) * other[1];
    }

  private:
    double phi_;
    const double spin_value_;
};

class HeisenbergSpin : public BaseSpin {
  public:
    explicit HeisenbergSpin(double spin_value) : spin_value_(spin_value) {}

    void randomize() override {
        phi_ = rng().uniform_real(0.0, 2.0 * M_PI);
        theta_ = rng().uniform_real(0.0, M_PI);
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

        return {spin_value_ * sin_theta * std::cos(phi_),
                spin_value_ * sin_theta * std::sin(phi_),
                spin_value_ * std::cos(theta_)};
    }
    double dot_product(const BaseSpin &other) const override {
        const auto &other_heisenberg =
            dynamic_cast<const HeisenbergSpin &>(other);
        double cos_angle = std::cos(phi_ - other_heisenberg.phi_) *
                               std::sin(theta_) *
                               std::sin(other_heisenberg.theta_) +
                           std::cos(theta_) * std::cos(other_heisenberg.theta_);

        return spin_value_ * spin_value_ * cos_angle;
    }
    double dot_product(const std::array<double, 3> &other) const override {
        const auto components = get_components();
        return components[0] * other[0] + components[1] * other[1] +
               components[2] * other[2];
    }

  private:
    double phi_, theta_;
    const double spin_value_;
};

class SpinFactory {
  public:
    static std::unique_ptr<BaseSpin> create(SpinModel model,
                                            double spin_value) {
        switch (model) {
        case SpinModel::ISING:
            return std::make_unique<IsingSpin>(spin_value);
        case SpinModel::XY:
            return std::make_unique<XYSpin>(spin_value);
        case SpinModel::HEISENBERG:
            return std::make_unique<HeisenbergSpin>(spin_value);
        default:
            throw std::invalid_argument("Unknown spin model");
        }
    }
};

} // namespace lattice

#endif // SPIN_HPP