#ifndef CP_ALGO_UTIL_CHECKPOINT_HPP
#define CP_ALGO_UTIL_CHECKPOINT_HPP
#include <iostream>
#include <chrono>
#include <string>
namespace cp_algo {
    template<bool final = false>
    void checkpoint([[maybe_unused]] std::string const& msg = "") {
#ifdef CP_ALGO_CHECKPOINT
        static double last = 0;
        double now = (double)clock() / CLOCKS_PER_SEC;
        double delta = now - last;
        last = now;
        if(msg.size()) {
            std::cerr << msg << ": " << (final ? now : delta) * 1000 << " ms\n";
        }
#endif
    }
}
#endif // CP_ALGO_UTIL_CHECKPOINT_HPP
