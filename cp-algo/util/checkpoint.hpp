#ifndef CP_ALGO_UTIL_CHECKPOINT_HPP
#define CP_ALGO_UTIL_CHECKPOINT_HPP
#include <iostream>
#include <chrono>
#include <string>
#include <map>
namespace cp_algo {
    template<bool final = false>
    void checkpoint([[maybe_unused]] auto const& msg = "") {
#ifdef CP_ALGO_CHECKPOINT
        static std::map<std::string, double> checkpoints;
        static double last = 0;
        double now = (double)clock() / CLOCKS_PER_SEC;
        double delta = now - last;
        last = now;
        if(msg.size() && !final) {
            checkpoints[msg] += delta;
        }
        if(final) {
            for(auto const& [key, value] : checkpoints) {
                std::cerr << key << ": " << value * 1000 << " ms\n";
            }
            std::cerr << "Total: " << now * 1000 << " ms\n";
        }
#endif
    }
}
#endif // CP_ALGO_UTIL_CHECKPOINT_HPP
