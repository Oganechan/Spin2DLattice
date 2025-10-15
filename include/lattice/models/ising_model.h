#pragma once

#include "base_model.h"

namespace lattice
{

    class IsingModel : public SpinModelBase
    {
    public:
        explicit IsingModel(const Geometry &geometry)
            : geometry_(geometry),
              num_atoms_(geometry.get_num_atoms())
        {
            magnetic_mask_.resize(num_atoms_, true);
            spins_.resize(num_atoms_);
        }

        void set_magnetic(int32_t atom_index, bool magnetic) override
        {
            if (atom_index < 0 || atom_index >= num_atoms_)
                throw std::invalid_argument("Index is out of range num_atoms");

            magnetic_mask_[atom_index] = magnetic;
        }

        bool is_magnetic(int32_t atom_index) const override
        {
            if (atom_index < 0 || atom_index >= num_atoms_)
                throw std::invalid_argument("Index is out of range num_atoms");

            return magnetic_mask_[atom_index];
        }

        const std::vector<bool> &get_magnetic_mask() const override
        {
            return magnetic_mask_;
        }

        void set_random_defects(double concentration) override
        {
            if (concentration < 0 || concentration > 1)
                throw std::invalid_argument("Incorrect concentration value. Concentration must be in [0,1]");

            std::fill(magnetic_mask_.begin(), magnetic_mask_.end(), true);

            int32_t deffect_count = static_cast<int32_t>((1.0 - concentration) * static_cast<double>(num_atoms_));

            std::vector<int32_t> indices(num_atoms_);
            for (int32_t i = 0; i < num_atoms_; ++i)
                indices[i] = i;
            std::shuffle(indices.begin(), indices.end(), Random::get_rng());

            for (int32_t i = 0; i < deffect_count; ++i)
                magnetic_mask_[indices[i]] = false;
        }

        void random_initialize() override
        {
            spins_.clear();
            spins_.reserve(num_atoms_);

            std::bernoulli_distribution dist(0.5);
            for (int32_t i = 0; i < num_atoms_; ++i)
                spins_.push_back(dist(Random::get_rng()) ? 1 : -1);
        }

        void ferromagnetic_initialize() override
        {
            spins_.clear();
            spins_.reserve(num_atoms_);

            for (int32_t i = 0; i < num_atoms_; ++i)
                spins_.push_back(1);
        }

        void antiferromagnetic_initialize() override
        {
            spins_.clear();
            spins_.reserve(num_atoms_);

            for (int32_t i = 0; i < num_atoms_; ++i)
                spins_[i] = (i % 2 == 0) ? 1 : -1;
        }

        SpinModel get_type() const override
        {
            return SpinModel::ISING;
        }

    private:
        const Geometry &geometry_;
        const int32_t num_atoms_;

        std::vector<int32_t> spins_;
        std::vector<bool> magnetic_mask_;
    };

}