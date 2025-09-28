#include "../include/config.h"

void Config::load(const std::string &filepath)
{
    std::unique_lock lock(mutex_);
    json new_config;

    std::ifstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Config file not found: " + filepath);
    file >> new_config;

    config_ = std::move(new_config);
}