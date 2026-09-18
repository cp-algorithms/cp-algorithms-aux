#include "cp-algo/math/fft.hpp"
#include <random>
#include <iostream>
#include <ranges>
using namespace cp_algo;
using namespace cp_algo::math;

// Naive reference over the first k coefficients.
template<typename T>
big_vector<T> naive(auto const& a, auto const& b, size_t k) {
    big_vector<T> c(k);
    for(size_t i = 0; i < std::size(a); i++) {
        for(size_t j = 0; j < std::size(b) && i + j < k; j++) {
            c[i + j] += a[i] * b[j];
        }
    }
    return c;
}

template<typename T>
big_vector<T> random_poly(size_t n, auto& rng) {
    big_vector<T> a(n);
    for(auto &x: a) {x = T(rng() % T::mod());}
    return a;
}

// One reusable transform serves several products, and every copy of one is a clone.
template<typename T>
void check_reuse(size_t as, size_t bs, auto& rng) {
    size_t need = as + bs - 1, cap = std::bit_ceil(need);
    auto a = random_poly<T>(as, rng), b = random_poly<T>(bs, rng), c = random_poly<T>(as, rng);
    auto ab = naive<T>(a, b, need), cb = naive<T>(c, b, need);

    auto B = fft::spectrum<T>(b, cap);
    big_vector<T> got(need);
    fft::spectrum<T>(a, cap).multiply(B, got, need);
    assert(got == ab);
    // The right operand survives, so the same transform multiplies another polynomial.
    got.assign(need, T(0));
    fft::spectrum<T>(c, cap).multiply(B, got, need);
    assert(got == cb);
    // A clone stands in for the consumed left operand.
    auto A = fft::spectrum<T>(a, cap);
    got.assign(need, T(0));
    A.clone().multiply(B, got, need);
    assert(got == ab);
    got.assign(need, T(0));
    std::move(A).multiply(B, got, need);
    assert(got == ab);
    // Truncated readback, and a product that fills the whole capacity.
    for(size_t k: {size_t(1), need / 2, need}) {
        if(!k) {continue;}
        big_vector<T> part(k);
        fft::spectrum<T>(a, cap).multiply(B, part, k);
        assert(std::equal(part.begin(), part.end(), ab.begin()));
    }
}

// Squares, views as input, and operands that use the wrapped half of a branch.
template<typename T>
void check_shapes(size_t n, auto& rng) {
    auto a = random_poly<T>(n, rng);
    size_t need = 2 * n - 1, cap = std::bit_ceil(need);
    auto aa = naive<T>(a, a, need);
    big_vector<T> got(need);
    fft::spectrum<T>(a, cap).square(got, need);
    assert(got == aa);
    // A view of the same coefficients builds the same transform.
    got.assign(need, T(0));
    auto view = a | std::views::take(n);
    fft::spectrum<T>(view, cap).multiply(fft::spectrum<T>(a, cap), got, need);
    assert(got == aa);
    // Operands longer than half the capacity exist only when the product still fits.
    size_t bs = std::max<size_t>(1, cap / 2 - n + 1);
    if(bs > 1 && n + bs - 1 <= cap) {
        auto b = random_poly<T>(bs, rng);
        auto want = naive<T>(a, b, n + bs - 1);
        big_vector<T> out(n + bs - 1);
        fft::spectrum<T>(a, cap).multiply(fft::spectrum<T>(b, cap), out, out.size());
        assert(out == want);
    }
}

// An accumulated sum of products needs one inverse transform for the whole sum.
template<typename T>
void check_product(size_t n, size_t terms, auto& rng) {
    size_t need = 2 * n - 1, cap = std::bit_ceil(need);
    big_vector<T> want(need);
    fft::product<T> acc(cap);
    big_vector<typename fft::product<T>::operand> xs, ys;
    for(size_t t = 0; t < terms; t++) {
        auto x = random_poly<T>(n, rng), y = random_poly<T>(n, rng);
        auto xy = naive<T>(x, y, need);
        for(size_t i = 0; i < need; i++) {want[i] += xy[i];}
        xs.emplace_back(x, cap);
        ys.emplace_back(y, cap);
    }
    // Each operand is read by several terms, as the multivariate product does.
    for(size_t t = 0; t < terms; t++) {acc.add(xs[t], ys[t]);}
    big_vector<T> got(need);
    std::move(acc).recover(got, need);
    assert(got == want);
}

// The cyclic product wraps every coefficient back into k positions.
template<typename T>
void check_cyclic(size_t k, auto& rng) {
    auto a = random_poly<T>(k, rng), b = random_poly<T>(k, rng);
    auto full = naive<T>(a, b, 2 * k - 1);
    big_vector<T> want(k);
    for(size_t i = 0; i < full.size(); i++) {want[i % k] += full[i];}
    auto got = a;
    fft::cyclic_mul(got, b, k);
    assert(got == want);
    // The same operand on both sides is a cyclic square.
    auto sq = a;
    auto fullsq = naive<T>(a, a, 2 * k - 1);
    big_vector<T> wantsq(k);
    for(size_t i = 0; i < fullsq.size(); i++) {wantsq[i % k] += fullsq[i];}
    fft::cyclic_mul(sq, sq, k);
    assert(sq == wantsq);
}

template<typename T>
void check(auto& rng) {
    for(size_t k: {size_t(8), size_t(64), size_t(256), size_t(1024)}) {
        check_cyclic<T>(k, rng);
    }
    for(auto [as, bs]: {std::pair{size_t(1), size_t(1)}, {size_t(8), size_t(5)},
                        {size_t(64), size_t(64)}, {size_t(100), size_t(7)},
                        {size_t(513), size_t(300)}, {size_t(1000), size_t(1000)}}) {
        check_reuse<T>(as, bs, rng);
    }
    for(size_t n: {size_t(1), size_t(5), size_t(64), size_t(129), size_t(600)}) {
        check_shapes<T>(n, rng);
    }
    for(size_t n: {size_t(4), size_t(64), size_t(300)}) {
        check_product<T>(n, 3, rng);
    }
}

int main() {
    std::mt19937_64 rng(20260918);
    check<modint<998244353>>(rng);      // d = 1
    check<modint<1000000007>>(rng);     // d = 5
    check<modint<2147483647>>(rng);     // the largest prime the rule admits
    check<modint<65537>>(rng);
    check<modint<int64_t(998244353)>>(rng);   // 8-byte storage, as the linear algebra uses
    check<modint<13>>(rng);
    check<modint<3>>(rng);
    for(int p: {998244353, 1000000007, 65537}) {
        dynamic_modint<>::with_mod(p, [&] {check<dynamic_modint<>>(rng);});
    }
    std::cout << "Reusable ring transforms and accumulated products passed\n";
}
