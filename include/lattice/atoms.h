#pragma once

#include "./models/ising_model.h"
#include "./models/heisenberg_model.h"

namespace lattice
{

    class Atoms
    {
    public:
        explicit Atoms(const Config &config);

        SpinModel get_type() { return model_->get_type(); }

        void set_magnetic(int32_t atom_index, bool magnetic) { model_->set_magnetic(atom_index, magnetic); }
        bool is_magnetic(int32_t atom_index) const { return model_->is_magnetic(atom_index); }
        const std::vector<bool> &get_magnetic_mask() { return model_->get_magnetic_mask(); }
        void set_random_defects(double concentration) { model_->set_random_defects(concentration); }

        int32_t count_magnetic_atoms() const { return model_->count_magnetic_atoms(); }
        int32_t count_defects() const { return model_->count_defects(); }
        std::vector<int32_t> get_magnetic_atom_indices() const { return model_->get_magnetic_atom_indices(); }
        std::vector<int32_t> get_defect_indices() const { return model_->get_defect_indices(); }

        SpinVariant get_spin(int32_t atom_index) const { return model_->get_spin(atom_index); }
        void set_spin(int32_t atom_index, const SpinVariant &value) { model_->set_spin(atom_index, value); }
        void flip_spins(const std::vector<int32_t> &indices) { model_->flip_spins(indices); }

        void random_initialize() { model_->random_initialize(); }

        const Geometry &geometry() const { return *geometry_; }
        const SpinModelBase &model() const { return *model_; }

    private:
        const Config &config_;
        std::unique_ptr<Geometry> geometry_;
        std::unique_ptr<SpinModelBase> model_;

        static SpinModel parse_model_type(const std::string &model_type_str);
        std::unique_ptr<SpinModelBase> create_model(SpinModel model_type);

        void initialize_geometry();
        void initialize_model();
    };

} // namespace lattice