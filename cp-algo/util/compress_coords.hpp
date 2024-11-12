#ifndef CP_ALGO_UTIL_COMPRESS_COORDS_HPP
#define CP_ALGO_UTIL_COMPRESS_COORDS_HPP
#include <algorithm>
#include <vector>
namespace cp_algo {
    std::vector<int> compress_coords(auto &coords) {
        static_assert(std::is_pointer_v<std::ranges::range_value_t<decltype(coords)>>);
        std::vector<int> original;
        original.reserve(size(coords));
        std::ranges::sort(coords, {}, [](int* x) {return *x;});
        int idx = -1, prev = -1;
        for(auto x: coords) {
            if(*x != prev) {
                idx++;
                prev = *x;
                original.push_back(*x);
            }
            *x = idx;
        }
        return original;
    }
}
#endif // CP_ALGO_UTIL_COMPRESS_COORDS_HPP