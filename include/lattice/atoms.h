#pragma once

#include "./models/ising_model.h"
#include "./models/heisenberg_model.h"
#include "./models/xy_model.h"

namespace lattice
{

    class Atoms
    {
    public:
        explicit Atoms(const Config &config);

        void set_magnetic(int32_t atom_index, bool magnetic) { model_->set_magnetic(atom_index, magnetic); }
        bool is_magnetic(int32_t atom_index) { return model_->is_magnetic(atom_index); }
        const std::vector<bool> &get_magnetic_mask() { return model_->get_magnetic_mask(); }
        void set_random_defects(double concentration) { model_->set_random_defects(concentration); }
        void random_initialize() { model_->random_initialize(); }
        void ferromagnetic_initialize() { model_->ferromagnetic_initialize(); }
        void antiferromagnetic_initialize() { model_->antiferromagnetic_initialize(); }
        SpinModel get_type() { return model_->get_type(); }

        const Geometry &geometry() const { return *geometry_; }
        SpinModelBase &model() { return *model_; }
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

}