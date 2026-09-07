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

size_t paired_products = 0, paired_gauss = 0;

template<typename base, typename row = modint_vec<base>>
void check_pairs() {
    using M = matrix<base, row>;
    std::mt19937 rng(623);
    std::array<size_t, 3> shapes[] = {
        {0, 0, 0}, {3, 0, 0}, {2, 9, 0}, {1, 7, 9}, {2, 9, 1},
        {3, 7, 7}, {4, 8, 4}, {5, 9, 5}, {8, 16, 17}, {9, 31, 31},
        {32, 33, 35}, {33, 65, 129}
    };
    for(auto [n, m, k]: shapes) for(int type = 0; type < 3; type++) {
        M a(n, m), b(m, k), expected(n, b.m());
        for(auto &x: a.elements()) x = type == 1 ? base(-1) : base(rng());
        for(auto &x: b.elements()) x = type == 1 ? base(-1) : base(rng());
        if(type == 2) {
            for(auto &x: a.elements()) if(rng() % 3) x = 0;
            for(auto &x: b.elements()) if(rng() % 3) x = 0;
        }
        for(size_t i = 0; i < n; i++)
        for(size_t j = 0; j < m; j++)
        for(size_t t = 0; t < b.m(); t++) expected[i][t] += a[i][j] * b[j][t];
        assert(a * b == expected);
        paired_products++;
    }
    for(size_t n: {1, 2, 3, 31, 32, 33, 65})
    for(size_t m: {0, 1, 3, 5, 33, 66, 101}) {
        M a(n, m);
        row source(m);
        for(auto &x: source) x = rng();
        for(size_t i = 0; i < n; i++) {
            for(auto &x: a[i]) x = rng();
            // Leave different deferred-reduction counts in adjacent rows.
            for(size_t j = 0; j < i % 13; j++) a[i].add_scaled(source, base(-1));
            if(i % 3 == 0 && m) a[i].normalize(i % m);
        }
        check_gauss<normal>(a);
        check_gauss<reverse>(a);
        paired_gauss += 2;
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
    check_pairs<modint<998244353LL>>();
    check_pairs<modint<1000000007LL>>();
    check_pairs<modint<1073741789LL>>();
    check_pairs<modint<998244353LL>, vec<modint<998244353LL>>>();
    for(int64_t p: {998244353LL, 1000000007LL, 1073741789LL}) {
        dynamic_modint<int64_t>::with_mod(p, [] {
            check_pairs<dynamic_modint<int64_t>>();
        });
    }
    std::cout << paired_products << " paired products and " << paired_gauss << " mixed-state Gaussian comparisons passed\n";
    std::cout << "96 accumulation cases, 936 Gaussian comparisons and rectangular application passed\n";
}
