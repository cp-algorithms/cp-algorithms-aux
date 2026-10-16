#include "cp-algo/number_theory/discrete_log.hpp"
#include <cassert>
#include <set>
using namespace cp_algo::math;
int main() {
    for(int m = 1; m <= 35; m++)
        for(int a = 0; a < m; a++)
            for(int b = 0; b < m; b++) {
                std::vector<std::optional<int>> first(m);
                int value = a;
                for(int x = 0; !first[value]; x++) {
                    first[value] = x;
                    value = int(int64_t(value) * b % m);
                }
                for(int c = 0; c < m; c++)
                    assert(discrete_log(b, c, m, a) == first[c]);
            }
    assert(discrete_log<int64_t>(5000000000LL, 0, 10000000000LL) == 2);
    assert(discrete_log<int64_t>(5000000000LL, 5000000000LL, 10000000000LL) == 1);
    assert(!discrete_log<int64_t>(5000000000LL, 3, 10000000000LL));
}
