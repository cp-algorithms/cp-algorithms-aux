#ifndef CP_ALGO_UTIL_COMPRESS_COORDS_HPP
#define CP_ALGO_UTIL_COMPRESS_COORDS_HPP
#include "sort.hpp"
#include <vector>
namespace cp_algo {
    // coords is a range of reference_wrapper<T>
    auto compress_coords(auto &coords) {
        std::vector<int> original;
        original.reserve(size(coords));
        radix_sort(coords);
        int idx = -1, prev = -1;
        for(auto &x: coords) {
            if(x != prev) {
                idx++;
                prev = x;
                original.push_back(x);
            }
            x.get() = idx;
        }
        return original;
    }
}
#endif // CP_ALGO_UTIL_COMPRESS_COORDS_HPP
