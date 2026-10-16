#include "cp-algo/number_theory/two_squares.hpp"
#include <cassert>
#include <set>
using namespace cp_algo::math;
template <class I> void check() {
    for(I n = 0; n <= 500; n++) {
        std::set<std::pair<int64_t, int64_t>> got, want;
        for(auto z : two_squares_all(n)) {
            static_assert(std::is_signed_v<typename decltype(z)::value_type>);
            assert(z.real() >= 0 && z.imag() >= 0);
            assert(__int128(z.real()) * z.real() + __int128(z.imag()) * z.imag() == n);
            got.emplace(z.real(), z.imag());
        }
        for(int64_t a = 0; a * a <= n; a++)
            for(int64_t b = 0; a * a + b * b <= n; b++)
                if(a * a + b * b == n)
                    want.emplace(a, b);
        assert(got == want);
    }
}
int main() {
    check<int>();
    check<int64_t>();
    check<uint32_t>();
    check<uint64_t>();
    auto z = two_squares_prime_any(uint32_t(4294967197U));
    assert(int64_t(z.real()) * z.real() + int64_t(z.imag()) * z.imag() == 4294967197LL);
    for(auto p : two_squares_all(uint64_t(1) << 62))
        assert(__int128(p.real()) * p.real() + __int128(p.imag()) * p.imag() == (__int128(1) << 62));
}
