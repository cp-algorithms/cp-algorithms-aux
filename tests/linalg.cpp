#include "cp-algo/linalg/matrix.hpp"
#include <random>
#include <iostream>

using namespace cp_algo::math;
using namespace cp_algo::linalg;

template<typename base>
void check_accumulation() {
    std::mt19937 rng(47);
    for(size_t n: {0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 65}) {
        modint_vec<base> a(n), b(n);
        std::vector<uint64_t> expected(n);
        for(size_t i = 0; i < n; i++) {
            a[i] = base::mod() - 1;
            b[i] = i % 2 ? base::mod() - 1 : rng() % base::mod();
            expected[i] = a[i].getr();
        }
        for(size_t t = 0; t < 1000; t++) {
            base scale = t % 2 ? base::mod() - 1 : rng() % base::mod();
            a.add_scaled(b, scale);
            for(size_t i = 0; i < n; i++) {
                expected[i] = (expected[i] + __uint128_t(b[i].getr()) * scale.getr()) % base::mod();
            }
            // Mix partial and full normalization without resetting the update sequence.
            if(n && t % 7 == 0) {
                size_t i = t % n;
                assert(a.normalize(i).getr() == expected[i]);
            }
            if(t % 23 == 0) {
                a.normalize();
                for(size_t i = 0; i < n; i++) assert(a[i].getr() == expected[i]);
            }
        }
        a.normalize();
        for(size_t i = 0; i < n; i++) assert(a[i].getr() == expected[i]);
    }
}

template<gauss_mode mode, typename M>
void check_gauss(M a) {
    M b = a;
    for(size_t i = 0; i < a.n(); i++) a.template eliminate<mode>(i);
    a.normalize();
    b.template gauss<mode>();
    assert(a == b);
}

void check_blocks() {
    using base = modint<998244353LL>;
    using M = matrix<base>;
    std::mt19937 rng(981);
    for(size_t n: {0, 1, 3, 15, 16, 17, 31, 32, 33, 63, 64, 65, 97}) {
        for(size_t m: {0, 1, 5, 17, 33, 66, 101, 513, 1025}) {
            for(int type = 0; type < 4; type++) {
                M a(n, m);
                if(type == 0) {
                    for(auto &x: a.elements()) x = rng();
                } else if(type == 1) {
                    for(size_t i = 0; i < std::min(n, m); i++) a[i][m - 1 - i] = rng();
                } else if(type == 2) {
                    M low(3, m);
                    for(auto &x: low.elements()) x = rng();
                    for(auto &row: a) for(auto &b: low) row.add_scaled(b, base(rng()));
                    a.normalize();
                } else {
                    for(auto &x: a.elements()) if(rng() % 20 == 0) x = rng();
                }
                check_gauss<normal>(a);
                check_gauss<reverse>(a);
            }
        }
    }
}

int main() {
    check_accumulation<modint<998244353LL>>();
    check_accumulation<modint<1000000007LL>>();
    check_accumulation<modint<1073741789LL>>();
    for(int64_t p: {998244353LL, 1000000007LL, 1073741789LL}) {
        dynamic_modint<int64_t>::with_mod(p, [] {
            check_accumulation<dynamic_modint<int64_t>>();
        });
    }
    matrix<int64_t, vec<int64_t>> a(3, 5);
    vec<int64_t> x(size_t(3));
    for(size_t i = 0; i < a.n(); i++) {
        x[i] = i + 1;
        for(size_t j = 0; j < a.m(); j++) a[i][j] = 5 * i + j;
    }
    auto y = a.apply(x);
    assert(y.size() == 5);
    for(size_t j = 0; j < y.size(); j++) assert(y[j] == 40 + 6 * int64_t(j));
    check_blocks();
    std::cout << "96 accumulation cases, 936 Gaussian comparisons and rectangular application passed\n";
}
