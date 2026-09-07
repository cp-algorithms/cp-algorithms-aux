#include "cp-algo/math/multivar_inv.hpp"
#include <random>
#include <iostream>
using namespace cp_algo::math;
template<typename T> void check() {
    std::mt19937 rng(917);
    for(size_t axes: {size_t(1), size_t(2), size_t(3)}) {
        for(size_t axis = 0; axis < axes; axis++) {
            for(size_t n: {size_t(1), size_t(2), size_t(31), size_t(65), size_t(129)}) {
                cp_algo::big_vector<size_t> dims(axes, 3), target(axes, 1);
                dims[axis] = n + 2; target[axis] = n;
                fft::multivar<T> a(dims);
                for(auto &x: a.data) {x = rng() % T::mod();}
                a.data[0] = 17;
                auto result = fft::multivar_inv(a, target);
                size_t stride = 1;
                for(size_t j = 0; j < axis; j++) {stride *= dims[j];}
                std::vector<T> want(n);
                want[0] = T(1) / a.data[0];
                for(size_t i = 1; i < n; i++) {
                    for(size_t j = 1; j <= i; j++) {want[i] -= a.data[stride * j] * want[i - j];}
                    want[i] *= want[0];
                }
                assert(result.dim == target && result.N == n);
                assert(std::ranges::equal(result.data, want));
            }
        }
    }
    fft::multivar<T> scalar(std::vector<size_t>{});
    scalar.data[0] = 19;
    auto result = fft::multivar_inv(scalar, std::vector<size_t>{});
    assert(result.N == 1 && result.data[0] * scalar.data[0] == T(1));
}
int main() {
    check<modint<998244353>>();
    check<modint<1000000007>>();
    std::cout << "Degenerate multivariate inverse properties passed under two primes\n";
}
