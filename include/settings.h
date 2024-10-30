#pragma once
#include <string>
#include <vector>

namespace settings {

    namespace random {
        inline auto seed         = 0ul;
    }

    namespace debug {
        inline bool allow_invalid_setting = true; //If this is true, then an invalid setting is simply ignored, if it is false then the program crashes.
    }

    namespace constants {
        //System sizes
        inline std::vector<int> sizes = {5, 5};

        //Number of simulation/thermalization steps
        inline long long unsigned int n_steps = 10;
        inline long long unsigned int n_therm = 0;

        inline double K = 1.0;
        inline std::vector<double> Ks = {0.0, 0.5, 1.0};
    }

    namespace io {
        inline std::string outfile = "simulation_data.h5";
        inline std::string settings_path = "../../settings/settings.json";
        //Probably want to set it to false in the final code.
        inline bool replace_file = true;
    }

    namespace log {
        inline long level = 2;
    }
}