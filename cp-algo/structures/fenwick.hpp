#ifndef CP_ALGO_STRUCTURES_FENWICK_HPP
#define CP_ALGO_STRUCTURES_FENWICK_HPP
#include <cassert>
#include <vector>
#include <bit>
namespace cp_algo::structures {
    template<typename T, typename Container = std::vector<T>>
    struct fenwick {
        size_t n;
        Container data;

        fenwick(auto &&range) {
            assign(range);
        }
        void to_prefix_sums() {
            for(size_t i = 1; i < n; i++) {
                if(i + (i & -i) <= n) {
                    data[i + (i & -i)] += data[i];
                }
            }
        }
        void assign(auto &&range) {
            n = size(range) - 1;
            data = move(range);
            to_prefix_sums();
        }
        void add(size_t x, T const& v) {
            for(++x; x <= n; x += x & -x) {
                data[x] += v;
            }
        }
        // sum of [0, r)
        T prefix_sum(size_t r) const {
            assert(r <= n);
            T res = 0;
            for(; r; r -= r & -r) {
                res += data[r];
            }
            return res;
        }
        // sum of [l, r)
        T range_sum(size_t l, size_t r) const {
            return prefix_sum(r) - prefix_sum(l);
        }
        // First r s.t. prefix_sum(r) >= k
        // Assumes data[x] >= 0 for all x
        size_t prefix_lower_bound(T k) const {
            int x = 0;
            for(size_t i = std::bit_floor(n); i; i /= 2) {
                if(x + i <= n && data[x + i] < k) {
                    k -= data[x + i];
                    x += i;
                }
            }
            return x;
        }
    };
}
#endif // CP_ALGO_STRUCTURES_FENWICK_HPP
