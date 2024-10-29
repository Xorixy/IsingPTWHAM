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
        inline int size_x = 4;
        inline int size_y = 4;
        inline int size_z = 1;

        inline std::vector<int> sizes;

        //Number of simulation/thermalization steps
        inline long long unsigned int n_steps = 10;
        inline long long unsigned int n_therm = 0;

        inline double K = 1.0;
    }

    namespace io {
        inline std::string outfile = "simulation.h5";
        inline std::string settings_path = "../../settings/settings.json";
        //Probably want to set it to false in the final code.
        inline bool replace_file = true;
    }

    namespace log {
        inline long level = 2;
    }
}