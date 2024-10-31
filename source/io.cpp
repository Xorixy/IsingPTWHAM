#include "../include/io.h"
#include <h5pp/details/h5ppEnums.h>
#include <spdlog/spdlog.h>

using json = nlohmann::json;

template<typename T>
std::vector<T> io::parse_array(const nlohmann::json &j, const std::string &name) {
    std::vector<T> vec;
    if (j.contains(name)) {
        try {
            for (auto& elem : j[name])
                vec.push_back(elem);
        } catch (std::exception &e) {
            if (!settings::debug::allow_invalid_setting) {
                const auto err_string = fmt::format("Invalid setting: {}. Error: {}", name, e.what());
                throw std::runtime_error(err_string);
            }
        }
    }
    return vec;
}

template<typename T>
std::optional<T> io::parse_setting(const nlohmann::json &j, const std::string &name) {
    if (j.contains(name)) {
        try {
            return j[name].get<T>();
        } catch (std::exception &e) {
            if (!settings::debug::allow_invalid_setting) {
                const auto err_string = fmt::format("Invalid setting: {}. Error: {}", name, e.what());
                throw std::runtime_error(err_string);
            }
        }
    }
    return std::nullopt;
}


void io::load_settings(const std::string& path) {
    std::ifstream f(path);
    json j_settings = json::parse(f);
    if (j_settings.contains("debug")) {
        auto j_debug = j_settings["debug"];
        if (std::optional<bool> set = parse_setting<bool>(j_debug, "allow_invalid_setting"); set.has_value()) {
            settings::debug::allow_invalid_setting = set.value();
        }
    }

    if (j_settings.contains("constants")) {
        auto j_constant = j_settings["constants"];
        if (std::vector<int> set = parse_array<int>(j_constant, "sizes"); !set.empty()) {
            settings::constants::sizes = set;
        }
        if (std::optional<long long unsigned int> set = parse_setting<long long unsigned int>(j_constant, "n_steps"); set.has_value()) {
            settings::constants::n_steps = set.value();
        }
        if (std::optional<long long unsigned int> set = parse_setting<long long unsigned int>(j_constant, "n_therm"); set.has_value()) {
            settings::constants::n_therm = set.value();
        }
        if (std::optional<long long unsigned int> set = parse_setting<long long unsigned int>(j_constant, "n_save"); set.has_value()) {
            settings::constants::n_save = set.value();
        }

        if (std::optional<double> set = parse_setting<double>(j_constant, "p_swap"); set.has_value()) {
            settings::constants::p_swap = set.value();
        }
        if (std::vector<double> set = parse_array<double>(j_constant, "Ks"); !set.empty()) {
            settings::constants::Ks = set;
        }
    }
    if (j_settings.contains("io")) {
        auto j_log = j_settings["io"];
        if (std::optional<std::string> set = parse_setting<std::string>(j_log, "outfile"); set.has_value()) {
            settings::io::outfile = set.value();
        }
    }
    if (j_settings.contains("log")) {
        auto j_log = j_settings["log"];
        if (std::optional<long> set = parse_setting<long>(j_log, "level"); set.has_value()) {
            settings::log::level = set.value();
        }
    }
}

h5pp::File io::try_to_open_file(const std::string& filename, bool readonly) {
    int n_tries = 0;
    h5pp::File file;
    while(true) {
        try {
            file = h5pp::File(filename, h5pp::FileAccess::READWRITE);
            break;
        } catch (std::exception &e) {
            fmt::print("Failed to open file with error: {}\n", e.what());
            fmt::print("Waiting 3 seconds and trying again...\n");
            std::this_thread::sleep_for(std::chrono::seconds(3));
            n_tries++;
            if (n_tries >= 10) throw e;
        }
    }
    return file;
}

/*
void io::save_settings(const std::string &path) {
    const json settings = {
        {"debug" ,
          {
              {"allow_invalid_setting" , settings::debug::allow_invalid_setting}
            }},
          {"worm" ,
            {
                {"p_single"   , settings::worm::p_single},
                {"p_double"   , settings::worm::p_double},
                {"p_quadruple", settings::worm::p_quadruple},

                {"p_free_move_single"   , settings::worm::p_free_move_single},
                {"p_free_kink_single"   , settings::worm::p_free_kink_single},
                {"p_snug_move_single"   , settings::worm::p_snug_move_single},
                {"p_snug_kink_single"   , settings::worm::p_snug_kink_single},
                {"p_snug_follow_single" , settings::worm::p_snug_follow_single},

                {"p_free_move_double"   , settings::worm::p_free_move_double},
                {"p_free_kink_double"   , settings::worm::p_free_kink_double},
                {"p_snug_move_double"   , settings::worm::p_snug_move_double},
                {"p_snug_kink_double"   , settings::worm::p_snug_kink_double},
                {"p_snug_follow_double" , settings::worm::p_snug_follow_double},

                {"p_free_move_quadruple"  , settings::worm::p_free_move_quadruple},
                {"p_free_kink_quadruple"  , settings::worm::p_free_kink_quadruple},
                {"p_snug_move_quadruple"  , settings::worm::p_snug_move_quadruple},
                {"p_snug_kink_quadruple"  , settings::worm::p_snug_kink_quadruple},
                {"p_snug_follow_quadruple", settings::worm::p_snug_follow_quadruple},

                {"C_single_snug" , settings::worm::C_single_snug},
                {"C_double_snug" , settings::worm::C_single_snug},
                {"C_quadruple_snug" , settings::worm::C_single_snug},
                {"C_single_free" , settings::worm::C_single_snug},
                {"C_double_free" , settings::worm::C_single_snug},
                {"C_quadruple_free" , settings::worm::C_single_snug}
              }},

            {"constants" ,
              {
                  {"size_x", settings::constants::size_x},
                  {"size_y", settings::constants::size_y},

                  {"n_steps", settings::constants::n_steps},
                  {"n_therm", settings::constants::n_therm},

                  {"beta", settings::constants::beta},
                  {"t"   , settings::constants::t},
                  {"d"   , settings::constants::d},
                  {"q"   , settings::constants::q},

                  {"mu", settings::constants::mu},
                  {"U" , settings::constants::U},
                  {"V" , settings::constants::V},

                  {"l_single"   , settings::constants::l_single},
                  {"l_double"   , settings::constants::l_double},
                  {"l_quadruple", settings::constants::l_quadruple},

                  {"fermionic_sign", settings::constants::fermionic_sign}
                }},
              {"save" ,
                {
                    {"windings", settings::save::windings},

                    {"correlations", settings::save::correlations},
                    {"annulus_size", settings::save::annulus_size},

                    {"time_series", settings::save::time_series},
                    {"save_interval", settings::save::save_interval}
                  }},
            {"log" ,
                {
                    {"level", settings::log::level},
                  }},
    };
    std::ofstream o(path);
    o << std::setw(4) << settings << std::endl;
}
*/
