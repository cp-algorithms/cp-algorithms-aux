#ifndef CP_ALGO_GEOMETRY_CLOSEST_PAIR_HPP
#define CP_ALGO_GEOMETRY_CLOSEST_PAIR_HPP
#include "../random/rng.hpp"
#include "../util/big_alloc.hpp"
#include "point.hpp"
#include <vector>
#include <map>
namespace cp_algo::geometry {
    // Rabin & Lipton
    auto closest_pair(auto const& r) {
        using point = std::decay_t<decltype(r[0])>;
        size_t n = size(r);
        int64_t md = 1e18;
        for(size_t i = 0; i < n / 100; i++) {
            auto A = random::rng() % n;
            auto B = random::rng() % n;
            if(A != B) {
                md = std::min(md, norm(r[A] - r[B]));
                if(md == 0) {
                    return std::pair{A, B};
                }
            }
        }
        std::map<point, big_vector<size_t>> neigs;
        md = (int64_t)ceil(sqrt((double)md));
        for(size_t i = 0; i < n; i++) {
            neigs[r[i] / md].push_back(i);
        }
        size_t a = 0, b = 1;
        md = norm(r[a] - r[b]);
        for(auto &[p, id]: neigs) {
            for(int dx: {-1, 0, 1}) {
                for(int dy: {-1, 0, 1}) {
                    auto pp = p + point{dx, dy};
                    if(!neigs.count(pp)) {
                        continue;
                    }
                    for(size_t i: neigs[pp]) {
                        for(size_t j: id) {
                            if(j == i) {
                                break;
                            }
                            int64_t cur = norm(r[i] - r[j]);
                            if(cur < md) {
                                md = cur;
                                a = i;
                                b = j;
                            }
                        }
                    }
                }
            }
        }
        return std::pair{a, b};
    }
}
#endif // CP_ALGO_GEOMETRY_CLOSEST_PAIR_HPP
