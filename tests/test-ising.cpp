#include "../include/ising.h"

int main(){
    std::vector<int> sizes;
    sizes.push_back(3);
    ising::Ising is(sizes, 1.0, 1);
    is.print_state();
    for (int in = 0 ; in < 10 ; in++) {
        is.run_sim_step(false);
        is.print_state();
    }
}