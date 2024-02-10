---
data:
  _extendedDependsOn:
  - icon: ':x:'
    path: cp-algo/algebra/common.hpp
    title: cp-algo/algebra/common.hpp
  _extendedRequiredBy:
  - icon: ':x:'
    path: cp-algo/algebra/fft.hpp
    title: cp-algo/algebra/fft.hpp
  - icon: ':warning:'
    path: cp-algo/algebra/matrix.hpp
    title: cp-algo/algebra/matrix.hpp
  - icon: ':x:'
    path: cp-algo/algebra/polynomial.hpp
    title: cp-algo/algebra/polynomial.hpp
  _extendedVerifiedWith:
  - icon: ':x:'
    path: verify/algebra/convolution107.test.cpp
    title: verify/algebra/convolution107.test.cpp
  _isVerificationFailed: true
  _pathExtension: hpp
  _verificationStatusIcon: ':x:'
  attributes:
    links:
    - https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm
  bundledCode: "#line 1 \"cp-algo/algebra/modular.hpp\"\n\n\n#line 1 \"cp-algo/algebra/common.hpp\"\
    \n\n\n#include <chrono>\n#include <random>\nnamespace algebra {\n    const int\
    \ maxn = 1 << 20;\n    const int magic = 250; // threshold for sizes to run the\
    \ naive algo\n    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());\
    \ \n\n    auto bpow(auto x, int64_t n, auto ans) {\n        for(; n; n /= 2, x\
    \ = x * x) {\n            if(n % 2) {\n                ans = ans * x;\n      \
    \      }\n        }\n        return ans;\n    }\n    template<typename T>\n  \
    \  T bpow(T const& x, int64_t n) {\n        return bpow(x, n, T(1));\n    }\n\n\
    \    template<typename T>\n    T fact(int n) {\n        static T F[maxn];\n  \
    \      static bool init = false;\n        if(!init) {\n            F[0] = T(1);\n\
    \            for(int i = 1; i < maxn; i++) {\n                F[i] = F[i - 1]\
    \ * T(i);\n            }\n            init = true;\n        }\n        return\
    \ F[n];\n    }\n    \n    template<typename T>\n    T rfact(int n) {\n       \
    \ static T F[maxn];\n        static bool init = false;\n        if(!init) {\n\
    \            F[maxn - 1] = T(1) / fact<T>(maxn - 1);\n            for(int i =\
    \ maxn - 2; i >= 0; i--) {\n                F[i] = F[i + 1] * T(i + 1);\n    \
    \        }\n            init = true;\n        }\n        return F[n];\n    }\n\
    \n    template<typename T>\n    T small_inv(int n) {\n        static T F[maxn];\n\
    \        static bool init = false;\n        if(!init) {\n            for(int i\
    \ = 1; i < maxn; i++) {\n                F[i] = rfact<T>(i) * fact<T>(i - 1);\n\
    \            }\n            init = true;\n        }\n        return F[n];\n  \
    \  }\n}\n\n#line 4 \"cp-algo/algebra/modular.hpp\"\n#include <algorithm>\n#include\
    \ <iostream>\n#include <optional>\nnamespace algebra {\n    template<int m>\n\
    \    struct modular {\n        // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm\n\
    \        // solves x^2 = y (mod m) assuming m is prime in O(log m).\n        //\
    \ returns std::nullopt if no sol.\n        std::optional<modular> sqrt() const\
    \ {\n            static modular y;\n            y = *this;\n            if(r ==\
    \ 0) {\n                return 0;\n            } else if(bpow(y, (m - 1) / 2)\
    \ != modular(1)) {\n                return std::nullopt;\n            } else {\n\
    \                while(true) {\n                    modular z = rng();\n     \
    \               if(z * z == *this) {\n                        return z;\n    \
    \                }\n                    struct lin {\n                       \
    \ modular a, b;\n                        lin(modular a, modular b): a(a), b(b)\
    \ {}\n                        lin(modular a): a(a), b(0) {}\n                \
    \        lin operator * (const lin& t) {\n                            return {\n\
    \                                a * t.a + b * t.b * y,\n                    \
    \            a * t.b + b * t.a\n                            };\n             \
    \           }\n                    } x(z, 1); // z + x\n                    x\
    \ = bpow(x, (m - 1) / 2);\n                    if(x.b != modular(0)) {\n     \
    \                   return x.b.inv();\n                    }\n               \
    \ }\n            }\n        }\n        \n        uint64_t r;\n        constexpr\
    \ modular(): r(0) {}\n        constexpr modular(int64_t rr): r(rr % m) {r = std::min<uint64_t>(r,\
    \ r + m);}\n        modular inv() const {return bpow(*this, m - 2);}\n       \
    \ modular operator - () const {return std::min(-r, m - r);}\n        modular operator\
    \ * (const modular &t) const {return r * t.r;}\n        modular operator / (const\
    \ modular &t) const {return *this * t.inv();}\n        modular& operator += (const\
    \ modular &t) {r += t.r; r = std::min<uint64_t>(r, r - m); return *this;}\n  \
    \      modular& operator -= (const modular &t) {r -= t.r; r = std::min<uint64_t>(r,\
    \ r + m); return *this;}\n        modular operator + (const modular &t) const\
    \ {return modular(*this) += t;}\n        modular operator - (const modular &t)\
    \ const {return modular(*this) -= t;}\n        modular& operator *= (const modular\
    \ &t) {return *this = *this * t;}\n        modular& operator /= (const modular\
    \ &t) {return *this = *this / t;}\n        \n        auto operator <=> (const\
    \ modular &t) const = default;\n        \n        explicit operator int() const\
    \ {return r;}\n        int64_t rem() const {return 2 * r > m ? r - m : r;}\n\n\
    \        static constexpr uint64_t mm = (uint64_t)m * m;\n        void add_unsafe(uint64_t\
    \ t) {r += t; r = std::min<uint64_t>(r, r - mm);}\n        modular& normalize()\
    \ {if(r >= m) r %= m; return *this;}\n    };\n    \n    template<int m>\n    std::istream&\
    \ operator >> (std::istream &in, modular<m> &x) {\n        return in >> x.r;\n\
    \    }\n    \n    template<int m>\n    std::ostream& operator << (std::ostream\
    \ &out, modular<m> const& x) {\n        return out << x.r % m;\n    }\n}\n\n"
  code: "#ifndef ALGEBRA_MODULAR_HPP\n#define ALGEBRA_MODULAR_HPP\n#include \"common.hpp\"\
    \n#include <algorithm>\n#include <iostream>\n#include <optional>\nnamespace algebra\
    \ {\n    template<int m>\n    struct modular {\n        // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm\n\
    \        // solves x^2 = y (mod m) assuming m is prime in O(log m).\n        //\
    \ returns std::nullopt if no sol.\n        std::optional<modular> sqrt() const\
    \ {\n            static modular y;\n            y = *this;\n            if(r ==\
    \ 0) {\n                return 0;\n            } else if(bpow(y, (m - 1) / 2)\
    \ != modular(1)) {\n                return std::nullopt;\n            } else {\n\
    \                while(true) {\n                    modular z = rng();\n     \
    \               if(z * z == *this) {\n                        return z;\n    \
    \                }\n                    struct lin {\n                       \
    \ modular a, b;\n                        lin(modular a, modular b): a(a), b(b)\
    \ {}\n                        lin(modular a): a(a), b(0) {}\n                \
    \        lin operator * (const lin& t) {\n                            return {\n\
    \                                a * t.a + b * t.b * y,\n                    \
    \            a * t.b + b * t.a\n                            };\n             \
    \           }\n                    } x(z, 1); // z + x\n                    x\
    \ = bpow(x, (m - 1) / 2);\n                    if(x.b != modular(0)) {\n     \
    \                   return x.b.inv();\n                    }\n               \
    \ }\n            }\n        }\n        \n        uint64_t r;\n        constexpr\
    \ modular(): r(0) {}\n        constexpr modular(int64_t rr): r(rr % m) {r = std::min<uint64_t>(r,\
    \ r + m);}\n        modular inv() const {return bpow(*this, m - 2);}\n       \
    \ modular operator - () const {return std::min(-r, m - r);}\n        modular operator\
    \ * (const modular &t) const {return r * t.r;}\n        modular operator / (const\
    \ modular &t) const {return *this * t.inv();}\n        modular& operator += (const\
    \ modular &t) {r += t.r; r = std::min<uint64_t>(r, r - m); return *this;}\n  \
    \      modular& operator -= (const modular &t) {r -= t.r; r = std::min<uint64_t>(r,\
    \ r + m); return *this;}\n        modular operator + (const modular &t) const\
    \ {return modular(*this) += t;}\n        modular operator - (const modular &t)\
    \ const {return modular(*this) -= t;}\n        modular& operator *= (const modular\
    \ &t) {return *this = *this * t;}\n        modular& operator /= (const modular\
    \ &t) {return *this = *this / t;}\n        \n        auto operator <=> (const\
    \ modular &t) const = default;\n        \n        explicit operator int() const\
    \ {return r;}\n        int64_t rem() const {return 2 * r > m ? r - m : r;}\n\n\
    \        static constexpr uint64_t mm = (uint64_t)m * m;\n        void add_unsafe(uint64_t\
    \ t) {r += t; r = std::min<uint64_t>(r, r - mm);}\n        modular& normalize()\
    \ {if(r >= m) r %= m; return *this;}\n    };\n    \n    template<int m>\n    std::istream&\
    \ operator >> (std::istream &in, modular<m> &x) {\n        return in >> x.r;\n\
    \    }\n    \n    template<int m>\n    std::ostream& operator << (std::ostream\
    \ &out, modular<m> const& x) {\n        return out << x.r % m;\n    }\n}\n#endif\
    \ // ALGEBRA_MODULAR_HPP\n"
  dependsOn:
  - cp-algo/algebra/common.hpp
  isVerificationFile: false
  path: cp-algo/algebra/modular.hpp
  requiredBy:
  - cp-algo/algebra/matrix.hpp
  - cp-algo/algebra/polynomial.hpp
  - cp-algo/algebra/fft.hpp
  timestamp: '2024-02-10 18:45:35+01:00'
  verificationStatus: LIBRARY_ALL_WA
  verifiedWith:
  - verify/algebra/convolution107.test.cpp
documentation_of: cp-algo/algebra/modular.hpp
layout: document
redirect_from:
- /library/cp-algo/algebra/modular.hpp
- /library/cp-algo/algebra/modular.hpp.html
title: cp-algo/algebra/modular.hpp
---
