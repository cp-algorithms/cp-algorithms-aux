#include "cp-algo/number_theory/two_squares.hpp"
#include <set>
using namespace cp_algo::math;

void check(int64_t n) {
    std::set<std::pair<int64_t, int64_t>> want, got;
    for(int64_t x = 0, y = int64_t(std::sqrt(n)); x * x <= n; x++) {
        while(x * x + y * y > n) {y--;}
        if(x * x + y * y == n) {want.emplace(x, y);}
    }
    for(auto z: two_squares_all(n)) {
        assert(z.real() >= 0 && z.imag() >= 0 && z.norm() == n);
        got.emplace(z.real(), z.imag());
    }
    assert(got == want);
}
int main() {
    cp_algo::random::gen.seed(71283);
    for(int64_t n = 0; n <= 2500; n++) {check(n);}
    for(int64_t n: {3125, 28125, 400000, 422500, 1267500}) {check(n);}
    std::cout << "Sums of two squares passed exhaustive enumeration checks\n";
}
