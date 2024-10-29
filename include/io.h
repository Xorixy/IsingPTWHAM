#pragma once
#include "settings.h"
#include <string>
#include <nlohmann/json.hpp>
#include <fstream>
#include <fmt/core.h>
#include <optional>

namespace io {
    void load_settings(const std::string& path);
    void save_settings(const std::string& path);
    template <typename T>
    std::optional<T> parse_setting(const nlohmann::json& j, const std::string& name);
}