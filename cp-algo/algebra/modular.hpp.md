---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/affine.hpp
    title: cp-algo/algebra/affine.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/common.hpp
    title: cp-algo/algebra/common.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/random/rng.hpp
    title: cp-algo/random/rng.hpp
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/fft.hpp
    title: cp-algo/algebra/fft.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/matrix.hpp
    title: cp-algo/algebra/matrix.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/polynomial.hpp
    title: cp-algo/algebra/polynomial.hpp
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/algebra/matrix/matrix_pow.test.cpp
    title: Pow of Matrix
  - icon: ':heavy_check_mark:'
    path: verify/algebra/polynomial/convolution107.test.cpp
    title: Convolution mod $10^9+7$
  - icon: ':heavy_check_mark:'
    path: verify/algebra/polynomial/find_linrec.test.cpp
    title: Find Linear Recurrence
  - icon: ':heavy_check_mark:'
    path: verify/algebra/polynomial/poly_sqrt.test.cpp
    title: Sqrt of Formal Power Series
  - icon: ':heavy_check_mark:'
    path: verify/algebra/polynomial/poly_sqrt.test.cpp
    title: Sqrt of Formal Power Series
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
    title: Dynamic Range Affine Range Sum
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links:
    - https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm
  bundledCode: "#line 1 \"cp-algo/algebra/modular.hpp\"\n\n\n#line 1 \"cp-algo/random/rng.hpp\"\
    \n\n\n#include <chrono>\n#include <random>\nnamespace cp_algo::random {\n    std::mt19937_64\
    \ rng(std::chrono::steady_clock::now().time_since_epoch().count()); \n}\n\n#line\
    \ 1 \"cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n#include <cassert>\n\
    namespace cp_algo::algebra {\n    template<typename base>\n    // a * x + b\n\
    \    struct lin {\n        base a = 1, b = 0;\n        std::optional<base> c;\n\
    \        lin() {}\n        lin(base b): a(0), b(b) {}\n        lin(base a, base\
    \ b): a(a), b(b) {}\n        lin(base a, base b, base _c): a(a), b(b), c(_c) {}\n\
    \n        // polynomial product modulo x^2 - c\n        lin operator * (const\
    \ lin& t) {\n            assert(c && t.c && *c == *t.c);\n            return {a\
    \ * t.b + b * t.a, b * t.b + a * t.a * (*c), *c};\n        }\n\n        // a *\
    \ (t.a * x + t.b) + b\n        lin apply(lin const& t) const {\n            return\
    \ {a * t.a, a * t.b + b};\n        }\n\n        void prepend(lin const& t) {\n\
    \            *this = t.apply(*this);\n        }\n\n        base eval(base x) const\
    \ {\n            return a * x + b;\n        }\n    };\n}\n\n#line 1 \"cp-algo/algebra/common.hpp\"\
    \n\n\n#include <cstdint>\nnamespace cp_algo::algebra {\n    const int maxn = 1\
    \ << 20;\n    const int magic = 250; // threshold for sizes to run the naive algo\n\
    \n    auto bpow(auto x, int64_t n, auto ans) {\n        for(; n; n /= 2, x = x\
    \ * x) {\n            if(n % 2) {\n                ans = ans * x;\n          \
    \  }\n        }\n        return ans;\n    }\n    template<typename T>\n    T bpow(T\
    \ const& x, int64_t n) {\n        return bpow(x, n, T(1));\n    }\n\n    template<typename\
    \ T>\n    T fact(int n) {\n        static T F[maxn];\n        static bool init\
    \ = false;\n        if(!init) {\n            F[0] = T(1);\n            for(int\
    \ i = 1; i < maxn; i++) {\n                F[i] = F[i - 1] * T(i);\n         \
    \   }\n            init = true;\n        }\n        return F[n];\n    }\n    \n\
    \    template<typename T>\n    T rfact(int n) {\n        static T F[maxn];\n \
    \       static bool init = false;\n        if(!init) {\n            F[maxn - 1]\
    \ = T(1) / fact<T>(maxn - 1);\n            for(int i = maxn - 2; i >= 0; i--)\
    \ {\n                F[i] = F[i + 1] * T(i + 1);\n            }\n            init\
    \ = true;\n        }\n        return F[n];\n    }\n\n    template<typename T>\n\
    \    T small_inv(int n) {\n        static T F[maxn];\n        static bool init\
    \ = false;\n        if(!init) {\n            for(int i = 1; i < maxn; i++) {\n\
    \                F[i] = rfact<T>(i) * fact<T>(i - 1);\n            }\n       \
    \     init = true;\n        }\n        return F[n];\n    }\n}\n\n#line 6 \"cp-algo/algebra/modular.hpp\"\
    \n#include <algorithm>\n#include <iostream>\n#line 9 \"cp-algo/algebra/modular.hpp\"\
    \nnamespace cp_algo::algebra {\n    template<int m>\n    struct modular {\n  \
    \      // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm\n        std::optional<modular>\
    \ sqrt() const {\n            if(r == 0) {\n                return 0;\n      \
    \      } else if(bpow(*this, (m - 1) / 2) != modular(1)) {\n                return\
    \ std::nullopt;\n            } else {\n                while(true) {\n       \
    \             modular z = random::rng();\n                    if(z * z == *this)\
    \ {\n                        return z;\n                    }\n              \
    \      lin<modular> x(1, z, *this); // x + z (mod x^2 - b)\n                 \
    \   x = bpow(x, (m - 1) / 2, lin<modular>(0, 1, *this));\n                   \
    \ if(x.a != modular(0)) {\n                        return x.a.inv();\n       \
    \             }\n                }\n            }\n        }\n        \n     \
    \   uint64_t r;\n        constexpr modular(): r(0) {}\n        constexpr modular(int64_t\
    \ rr): r(rr % m) {r = std::min<uint64_t>(r, r + m);}\n        modular inv() const\
    \ {return bpow(*this, m - 2);}\n        modular operator - () const {return std::min(-r,\
    \ m - r);}\n        modular operator * (const modular &t) const {return r * t.r;}\n\
    \        modular operator / (const modular &t) const {return *this * t.inv();}\n\
    \        modular& operator += (const modular &t) {r += t.r; r = std::min<uint64_t>(r,\
    \ r - m); return *this;}\n        modular& operator -= (const modular &t) {r -=\
    \ t.r; r = std::min<uint64_t>(r, r + m); return *this;}\n        modular operator\
    \ + (const modular &t) const {return modular(*this) += t;}\n        modular operator\
    \ - (const modular &t) const {return modular(*this) -= t;}\n        modular& operator\
    \ *= (const modular &t) {return *this = *this * t;}\n        modular& operator\
    \ /= (const modular &t) {return *this = *this / t;}\n        \n        auto operator\
    \ <=> (const modular &t) const = default;\n        \n        explicit operator\
    \ int() const {return r;}\n        int64_t rem() const {return 2 * r > m ? r -\
    \ m : r;}\n\n        static constexpr uint64_t mm = (uint64_t)m * m;\n       \
    \ void add_unsafe(uint64_t t) {r += t; r = std::min<uint64_t>(r, r - mm);}\n \
    \       modular& normalize() {if(r >= m) r %= m; return *this;}\n    };\n    \n\
    \    template<int m>\n    std::istream& operator >> (std::istream &in, modular<m>\
    \ &x) {\n        return in >> x.r;\n    }\n    \n    template<int m>\n    std::ostream&\
    \ operator << (std::ostream &out, modular<m> const& x) {\n        return out <<\
    \ x.r % m;\n    }\n}\n\n"
  code: "#ifndef CP_ALGO_ALGEBRA_MODULAR_HPP\n#define CP_ALGO_ALGEBRA_MODULAR_HPP\n\
    #include \"../random/rng.hpp\"\n#include \"affine.hpp\"\n#include \"common.hpp\"\
    \n#include <algorithm>\n#include <iostream>\n#include <optional>\nnamespace cp_algo::algebra\
    \ {\n    template<int m>\n    struct modular {\n        // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm\n\
    \        std::optional<modular> sqrt() const {\n            if(r == 0) {\n   \
    \             return 0;\n            } else if(bpow(*this, (m - 1) / 2) != modular(1))\
    \ {\n                return std::nullopt;\n            } else {\n            \
    \    while(true) {\n                    modular z = random::rng();\n         \
    \           if(z * z == *this) {\n                        return z;\n        \
    \            }\n                    lin<modular> x(1, z, *this); // x + z (mod\
    \ x^2 - b)\n                    x = bpow(x, (m - 1) / 2, lin<modular>(0, 1, *this));\n\
    \                    if(x.a != modular(0)) {\n                        return x.a.inv();\n\
    \                    }\n                }\n            }\n        }\n        \n\
    \        uint64_t r;\n        constexpr modular(): r(0) {}\n        constexpr\
    \ modular(int64_t rr): r(rr % m) {r = std::min<uint64_t>(r, r + m);}\n       \
    \ modular inv() const {return bpow(*this, m - 2);}\n        modular operator -\
    \ () const {return std::min(-r, m - r);}\n        modular operator * (const modular\
    \ &t) const {return r * t.r;}\n        modular operator / (const modular &t) const\
    \ {return *this * t.inv();}\n        modular& operator += (const modular &t) {r\
    \ += t.r; r = std::min<uint64_t>(r, r - m); return *this;}\n        modular& operator\
    \ -= (const modular &t) {r -= t.r; r = std::min<uint64_t>(r, r + m); return *this;}\n\
    \        modular operator + (const modular &t) const {return modular(*this) +=\
    \ t;}\n        modular operator - (const modular &t) const {return modular(*this)\
    \ -= t;}\n        modular& operator *= (const modular &t) {return *this = *this\
    \ * t;}\n        modular& operator /= (const modular &t) {return *this = *this\
    \ / t;}\n        \n        auto operator <=> (const modular &t) const = default;\n\
    \        \n        explicit operator int() const {return r;}\n        int64_t\
    \ rem() const {return 2 * r > m ? r - m : r;}\n\n        static constexpr uint64_t\
    \ mm = (uint64_t)m * m;\n        void add_unsafe(uint64_t t) {r += t; r = std::min<uint64_t>(r,\
    \ r - mm);}\n        modular& normalize() {if(r >= m) r %= m; return *this;}\n\
    \    };\n    \n    template<int m>\n    std::istream& operator >> (std::istream\
    \ &in, modular<m> &x) {\n        return in >> x.r;\n    }\n    \n    template<int\
    \ m>\n    std::ostream& operator << (std::ostream &out, modular<m> const& x) {\n\
    \        return out << x.r % m;\n    }\n}\n#endif // CP_ALGO_ALGEBRA_MODULAR_HPP\n"
  dependsOn:
  - cp-algo/random/rng.hpp
  - cp-algo/algebra/affine.hpp
  - cp-algo/algebra/common.hpp
  isVerificationFile: false
  path: cp-algo/algebra/modular.hpp
  requiredBy:
  - cp-algo/algebra/matrix.hpp
  - cp-algo/algebra/polynomial.hpp
  - cp-algo/algebra/fft.hpp
  timestamp: '2024-02-11 00:23:03+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  - verify/algebra/matrix/matrix_pow.test.cpp
  - verify/algebra/polynomial/convolution107.test.cpp
  - verify/algebra/polynomial/find_linrec.test.cpp
  - verify/algebra/polynomial/poly_sqrt.test.cpp
  - verify/algebra/polynomial/poly_sqrt.test.cpp
documentation_of: cp-algo/algebra/modular.hpp
layout: document
redirect_from:
- /library/cp-algo/algebra/modular.hpp
- /library/cp-algo/algebra/modular.hpp.html
title: cp-algo/algebra/modular.hpp
---
