#pragma once
#include "settings.h"
#include <string>
#include <nlohmann/json.hpp>
#include <fstream>
#include <fmt/core.h>
#include <optional>
#include <h5pp/h5pp.h>


namespace io {
    inline h5pp::File outfile;

    void load_settings(const std::string& path);
    void save_settings(const std::string& path);

    template <typename T>
    std::optional<T> parse_setting(const nlohmann::json& j, const std::string& name);

    template <typename T>
    std::vector<T> parse_array(const nlohmann::json& j, const std::string& name);

    template <typename T>
    void save_vector(std::vector<T> & vector, const std::string& name) {
        io::outfile.writeDataset<std::vector<T>>(vector, name);
    }

    h5pp::File try_to_open_file(const std::string& filename, bool readonly);
}