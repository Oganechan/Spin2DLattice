#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "../../external/nlohmann/json.hpp"
#include <fstream>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#include <string>

using json = nlohmann::json;

class Config {
  private:
    mutable std::shared_mutex mutex_;
    json config_ = json::object();

  public:
    void load(const std::string &filepath) {
        std::unique_lock lock(mutex_);
        json new_config;

        std::ifstream file(filepath);
        if (!file.is_open())
            throw std::runtime_error("Config file not found: " + filepath);
        file >> new_config;

        config_ = std::move(new_config);
    }

    template <typename T>
    T get(const std::string &key, const T &default_value = T()) const {
        std::shared_lock lock(mutex_);

        json current = config_;
        size_t pos = 0;
        std::string k = key;

        while ((pos = k.find('.')) != std::string::npos) {
            std::string part = k.substr(0, pos);

            if (!current.contains(part))
                return default_value;

            current = current[part];
            k = k.substr(pos + 1);
        }

        return current.value(k, default_value);
    }

    template <typename T> void set(const std::string &key, const T &value) {
        std::unique_lock lock(mutex_);

        json *current = &config_;
        size_t pos = 0;
        std::string k = key;

        while ((pos = k.find('.')) != std::string::npos) {
            std::string part = k.substr(0, pos);
            if (!current->contains(part))
                (*current)[part] = json::object();
            current = &(*current)[part];
            k = k.substr(pos + 1);
        }
        (*current)[k] = value;
    }
};

#endif // CONFIG_HPP