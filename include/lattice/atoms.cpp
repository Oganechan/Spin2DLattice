#include "atoms.h"

lattice::Atoms::Atoms(const Config &config)
    : config_(config)
{
    initialize_geometry();
    initialize_model();
}

lattice::SpinModel lattice::Atoms::parse_model_type(const std::string &model_type_str)
{
    static const std::unordered_map<std::string, SpinModel> mapping = {
        {"ising", SpinModel::ISING},
        {"heisenberg", SpinModel::HEISENBERG},
        {"xy", SpinModel::XY}};

    if (auto it = mapping.find(model_type_str); it != mapping.end())
        return it->second;

    throw std::invalid_argument("Unknown model type: " + model_type_str);
}

std::unique_ptr<lattice::SpinModelBase> lattice::Atoms::create_model(SpinModel model_type)
{
    switch (model_type)
    {
    case SpinModel::ISING:
        return std::make_unique<IsingModel>(*geometry_);

    case SpinModel::HEISENBERG:
        // return std::make_unique<HeisenbergModel>(*geometry_);
        throw;

    case SpinModel::XY:
        // return std::make_unique<XYModel>(*geometry_);
        throw;

    default:
        throw std::invalid_argument("Unknown spin model type");
    }
}

void lattice::Atoms::initialize_geometry()
{
    geometry_ = std::make_unique<Geometry>(config_);
}

void lattice::Atoms::initialize_model()
{
    std::string model_type_str = config_.get<std::string>("physical.spin_model");
    SpinModel model_type = parse_model_type(model_type_str);

    model_ = create_model(model_type);
}