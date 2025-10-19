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

        SpinModel get_type() const override { return SpinModel::ISING; }

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

        const std::vector<bool> &get_magnetic_mask() const override { return magnetic_mask_; }

        void set_random_defects(double concentration) override
        {
            if (concentration < 0 || concentration > 1)
                throw std::invalid_argument("Incorrect concentration value. Concentration must be in [0,1]");

            std::fill(magnetic_mask_.begin(), magnetic_mask_.end(), true);

            int32_t defect_count = static_cast<int32_t>(concentration * static_cast<double>(num_atoms_));

            std::vector<int32_t> indices(num_atoms_);
            for (int32_t i = 0; i < num_atoms_; ++i)
                indices[i] = i;
            std::shuffle(indices.begin(), indices.end(), Random::get_rng());

            for (int32_t i = 0; i < defect_count; ++i)
                magnetic_mask_[indices[i]] = false;
        }

        int32_t count_magnetic_atoms() const override
        {
            int32_t count = 0;
            for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
                count += magnetic_mask_[i];

            return count;
        }

        int32_t count_defects() const override
        {
            int32_t count = 0;
            for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
                count += !magnetic_mask_[i];

            return count;
        }

        std::vector<int32_t> get_magnetic_atom_indices() const override
        {
            std::vector<int32_t> indices;
            indices.reserve(magnetic_mask_.size());
            for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
                if (magnetic_mask_[i])
                    indices.push_back(i);

            return indices;
        }

        std::vector<int32_t> get_defect_indices() const override
        {
            std::vector<int32_t> indices;
            indices.reserve(magnetic_mask_.size());
            for (int32_t i = 0; i < magnetic_mask_.size(); ++i)
                if (!magnetic_mask_[i])
                    indices.push_back(i);

            return indices;
        }

        SpinVariant get_spin(int32_t atom_index) const override
        {
            if (atom_index < 0 || atom_index >= num_atoms_)
                throw std::invalid_argument("Index is out of range num_atoms");

            if (!magnetic_mask_[atom_index])
                throw std::invalid_argument("Atom is not magnetic");

            return spins_[atom_index];
        }

        void set_spin(int32_t atom_index, const SpinVariant &value) override
        {
            if (atom_index < 0 || atom_index >= num_atoms_)
                throw std::invalid_argument("Index is out of range num_atoms");

            if (!magnetic_mask_[atom_index])
                throw std::invalid_argument("Cannot set spin for non-magnetic atom");

            if (auto *spin_value = std::get_if<int32_t>(&value))
            {
                if (*spin_value != -1 && *spin_value != 1)
                    throw std::invalid_argument("Ising spin must be -1 or +1");
                spins_[atom_index] = *spin_value;
            }
            else
                throw std::invalid_argument("Invalid spin type for Ising model");
        }

        void flip_spins(const std::vector<int32_t> &indices) override
        {
            for (int32_t idx : indices)
            {
                if (idx < 0 || idx >= num_atoms_)
                    throw std::invalid_argument("Index is out of range num_atoms");

                if (!magnetic_mask_[idx])
                    throw std::invalid_argument("Cannot flip spin of non-magnetic atom");

                spins_[idx] *= -1;
            }
        }

        void random_initialize() override
        {
            spins_.clear();
            spins_.reserve(num_atoms_);

            std::bernoulli_distribution dist(0.5);
            for (int32_t i = 0; i < num_atoms_; ++i)
                spins_.push_back(dist(Random::get_rng()) ? 1 : -1);
        }

    private:
        const Geometry &geometry_;
        const int32_t num_atoms_;

        std::vector<int32_t> spins_;
        std::vector<bool> magnetic_mask_;
    };

} // namespace lattice