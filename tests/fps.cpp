#include "cp-algo/math/fps.hpp"
#include "cp-algo/math/poly.hpp"
#include <random>
#include <iostream>
using namespace cp_algo;
using namespace cp_algo::math;
template<typename T> void products() {
    using P = poly_t<T>;
    std::mt19937 rng(42);
    for(int m: {0, 1, 2, 63, 64, 65, 127, 128, 129, 257}) {
        typename P::Vector a(521), b(m), want(a.size() + b.size());
        for(auto &x: a) {x = rng() % 11;}
        for(auto &x: b) {x = rng() % 11;}
        for(size_t i = 0; i < a.size(); i++) {
            for(size_t j = 0; j < b.size(); j++) {want[i + j] += a[i] * b[j];}
        }
        size_t called = 0, allowed = 0;
        fps<T> dynamic([&](size_t n, auto const&) {
            assert(n == called++ && n <= allowed);
            return n < a.size() ? a[n] : T(0);
        });
        fps<T> fixed{P(b)};
        auto left = fixed * dynamic, right = dynamic * P(b);
        auto temporary = fps<T>(P(b)) * dynamic;
        assert(fixed[1000000] == T(0));
        assert(fixed.prefix(m + 100) == P(b));
        fixed = fps<T>();
        for(size_t n = 0; n < want.size() + 10; n++) {
            allowed = n;
            T expected = n < want.size() ? want[n] : T(0);
            assert(left[n] == expected && right[n] == expected && temporary[n] == expected);
        }
        assert(called <= want.size() + 10);
        assert(left.prefix(17) == P(typename P::Vector(want.begin(), want.begin() + 17)));
        fps<T> dynamic_b([&](size_t n, auto const&) {return n < b.size() ? b[n] : T(0);});
        auto full = dynamic * dynamic_b;
        assert(full.prefix(want.size()) == P(want));
    }
    // The changing input depends on the preceding output; lookahead would form a cycle.
    for(bool known: {false, true}) {
        P kernel({1, 2, 3});
        fps<T> product;
        fps<T> input([&](size_t n, auto const&) {return n ? product[n - 1] : T(1);});
        fps<T> factor = known ? fps<T>(kernel) : fps<T>([kernel](size_t n, auto const&) {
            return kernel[int(n)];
        });
        product = factor * input;
        typename P::Vector a, expected;
        for(size_t n = 0; n < 20; n++) {
            a.push_back(n ? expected.back() : T(1));
            T c = 0;
            for(size_t j = 0; j <= std::min(n, size_t(2)); j++) {c += kernel.a[j] * a[n - j];}
            expected.push_back(c);
            assert(product[n] == c);
        }
    }
}
template<typename T> void powers() {
    using P = poly_t<T>;
    std::mt19937 rng(713);
    for(int n: {0, 1, 2, 31, 32, 33, 63, 64, 65, 127, 128, 129, 257}) {
        for(int shift: {0, 1, 7, n}) {
            typename P::Vector a(n);
            for(int i = shift; i < n; i++) {a[i] = 1 + rng() % (T::mod() - 1);}
            for(int64_t k: {0LL, 1LL, 2LL, 7LL, 65LL, 1000000000000000000LL}) {
                P input(a), want = pow(input, k, n);
                auto fixed = pow(fps<T>(input), k);
                size_t allowed = 0, calls = 0;
                auto generated = fps<T>([&](size_t i, auto const&) {
                    assert(i == calls++ && i <= allowed);
                    return i < a.size() ? a[i] : T(0);
                });
                auto lazy = pow(generated, k);
                for(size_t i = 0; i < size_t(n); i++) {
                    allowed = i;
                    assert(fixed[i] == want[int(i)] && lazy[i] == want[int(i)]);
                }
                assert(fixed.prefix(n) == want && lazy.prefix(n) == want);
                if(k == 0) {assert(calls == 0);}
            }
        }
    }
    // Known sparse factors exercise both sides of the product dispatch cutoff.
    for(int terms: {1, 16, 17}) {
        typename P::Vector a(401);
        a[0] = 3;
        for(int j = 1; j <= terms; j++) {a[j * 19] = j + 7;}
        assert(pow(fps<T>(P(a)), 123456789).prefix(521) == pow(P(a), 123456789, 521));
    }
    bool threw = false;
    try {pow(fps<T>(T(1)), -1);} catch(std::domain_error const&) {threw = true;}
    assert(threw);
}
int main() {
    powers<modint<998244353>>();
    powers<modint<1000000007>>();
    products<modint<998244353>>();
    products<modint<1000000007>>();
    products<long long>();
    using T = modint<998244353>;
    using P = poly_t<T>;
    std::mt19937 rng(7);
    for(int terms: {0, 1, 15, 16, 17}) {
        P::Vector a(401);
        a[0] = 1;
        for(int j = 1; j <= terms; j++) {a[j * 19] = rng() % T::mod();}
        auto input = fps<T>(P(a));
        auto stream = fps<T>([](size_t n, auto const&) {return T(n + 1);});
        P::Vector b(521);
        for(size_t i = 0; i < b.size(); i++) {b[i] = T(i + 1);}
        assert((input * stream).prefix(521) == (P(a) * P(b)).mod_xk(521));
        assert(inv(input).prefix(521) == inv(P(a), 521));
        assert(log(input).prefix(521) == log(P(a), 521));
        a[0] = 0;
        assert(exp(fps<T>(P(a))).prefix(521) == exp(P(a), 521));
    }
    // The cached divisor table is finite; lazy evaluation may continue beyond it.
    fps<T> ones([](size_t, auto const&) {return T(1);});
    auto integral = integr(ones), logarithm = log(fps<T>(P({1, 1})));
    auto exponential = exp(fps<T>(P({0, 1})));
    T factorial = 1;
    for(size_t n = 1; n <= size_t(maxn) + 1; n++) {
        factorial *= T(n);
        if(n + 1 >= size_t(maxn)) {
            assert(integral[n] == T(1) / T(n));
            assert(logarithm[n] == T(n % 2 ? 1 : -1) / T(n));
            assert(exponential[n] * factorial == T(1));
        }
    }
    auto integer_integral = integr(fps<long long>([](size_t, auto const&) {return 6LL;}));
    assert(integer_integral[2] == 3 && integer_integral[4] == 1);
    auto check_integral = []<typename U>() {
        auto q = integr(fps<U>([](size_t, auto const&) {return U(1);}));
        assert(q[2] == U(1) / U(2));
    };
    check_integral.template operator()<modint<17>>();
    for(int mod: {101, 103}) {
        dynamic_modint<>::with_mod(mod, [&] {check_integral.template operator()<dynamic_modint<>>();});
    }
    for(int m: {1, 2, 63, 64, 65, 129, 257}) {
        P::Vector a(m);
        for(auto &x: a) {x = rng() % T::mod();}
        a[0] = 1;
        for(int kind = 0; kind < 3; kind++) {
            a[0] = kind == 2 ? 0 : 1;
            size_t allowed = 0, calls = 0;
            fps<T> changing([&](size_t n, auto const&) {
                assert(n == calls++ && n <= allowed);
                return n < a.size() ? a[n] : T(0);
            });
            fps<T> known{P(a)};
            auto apply = [kind](auto p) {return kind == 0 ? inv(p) : kind == 1 ? log(p) : exp(p);};
            auto semi = apply(known), full = apply(changing);
            P expected = kind == 0 ? inv(P(a), 521) : kind == 1 ? log(P(a), 521) : exp(P(a), 521);
            for(size_t n = 0; n < 521; n++) {
                allowed = n;
                assert(semi[n] == expected[int(n)] && full[n] == expected[int(n)]);
            }
            assert(calls == 521);
        }
    }
    fps<T> cycle;
    cycle = fps<T>([&](size_t n, auto const&) {return cycle[n + 1];});
    for(int attempt = 0; attempt < 2; attempt++) {
        bool rejected = false;
        try {(void)cycle[0];} catch(std::logic_error const&) {rejected = true;}
        assert(rejected);
    }
    int attempts = 0;
    fps<T> retry([&](size_t n, auto const& prev) {
        if(n == 3 && attempts++ == 0) {throw std::runtime_error("retry");}
        return prev.empty() ? T(1) : prev.back() + T(1);
    });
    bool rejected = false;
    try {(void)retry[5];} catch(std::runtime_error const&) {rejected = true;}
    assert(rejected && retry[2] == T(3) && retry[5] == T(6));
    std::cout << "FPS fixed, dynamic, cached, and causal products passed\n";
}
