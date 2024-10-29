#pragma once
#include <vector>
#include <fmt/core.h>


namespace stats {

    template <typename T>
    requires std::is_arithmetic_v<T>
    double relaxation_time(std::vector<T> & series, double t_frac = 0.1) {
        double mean = 0.0;
        double mean_sqr = 0.0;
        for (auto value : series) {
            mean += value;
            mean_sqr += value * value;
        }
        if (mean*mean == mean_sqr) {
            return 0.0;
        }
        mean /= series.size();
        mean_sqr /= series.size();
        double relaxation_time = 0.0;
        for (int t = 0 ; t < round(t_frac*(series.size() - 1)); t++) {
            double phi_t = 0.0;
            for (int in = 0 ; in < series.size() - t ; in++) {
                phi_t += series[in] * series[in + t];
            }
            phi_t /= series.size() - t;
            phi_t -= mean*mean;
            phi_t /= mean_sqr - mean*mean;
            fmt::print("{}\n", phi_t);
            //if (t % static_cast<int>(round(0.01*t_frac*(series.size() - 1))) == 0) fmt::print("{}\n", t);
            if (phi_t < 0) break;
            relaxation_time += phi_t;
        }
        return relaxation_time;
    }


}