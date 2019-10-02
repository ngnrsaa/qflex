#ifndef QFLEX_INPUT_
#define QFLEX_INPUT_

#include <iostream>
#include <string>

namespace qflex{
    struct QflexInput {
        int grid_height;
        int grid_width;
        int super_cycles;

        //deprecated?
        double fidelity;

        bool enable_timing_logs = false;

        std::istream* circuit_data;
        std::istream* ordering_data;
        std::istream* grid_data;
        std::string initial_state;
        std::string final_state_A;
    };
}

#endif