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
    path: cp-algo/algebra/matrix.hpp
    title: cp-algo/algebra/matrix.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/random/rng.hpp
    title: cp-algo/random/rng.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/pow_of_matrix
    document_title: Pow of Matrix
    links:
    - https://judge.yosupo.jp/problem/pow_of_matrix
  bundledCode: "#line 1 \"verify/algebra/matrix/matrix_pow.test.cpp\"\n// @brief Pow\
    \ of Matrix\n#define PROBLEM \"https://judge.yosupo.jp/problem/pow_of_matrix\"\
    \n#line 1 \"cp-algo/algebra/matrix.hpp\"\n\n\n#line 1 \"cp-algo/algebra/common.hpp\"\
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
    \     init = true;\n        }\n        return F[n];\n    }\n}\n\n#line 1 \"cp-algo/algebra/modular.hpp\"\
    \n\n\n#line 1 \"cp-algo/random/rng.hpp\"\n\n\n#include <chrono>\n#include <random>\n\
    namespace cp_algo::random {\n    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());\
    \ \n}\n\n#line 1 \"cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n#include\
    \ <cassert>\nnamespace cp_algo::algebra {\n    // a * x + b\n    template<typename\
    \ base>\n    struct lin {\n        base a = 1, b = 0;\n        std::optional<base>\
    \ c;\n        lin() {}\n        lin(base b): a(0), b(b) {}\n        lin(base a,\
    \ base b): a(a), b(b) {}\n        lin(base a, base b, base _c): a(a), b(b), c(_c)\
    \ {}\n\n        // polynomial product modulo x^2 - c\n        lin operator * (const\
    \ lin& t) {\n            assert(c && t.c && *c == *t.c);\n            return {a\
    \ * t.b + b * t.a, b * t.b + a * t.a * (*c), *c};\n        }\n\n        // a *\
    \ (t.a * x + t.b) + b\n        lin apply(lin const& t) const {\n            return\
    \ {a * t.a, a * t.b + b};\n        }\n\n        void prepend(lin const& t) {\n\
    \            *this = t.apply(*this);\n        }\n\n        base eval(base x) const\
    \ {\n            return a * x + b;\n        }\n    };\n\n    // (ax+b) / (cx+d)\n\
    \    template<typename base>\n    struct linfrac {\n        // default constructor\
    \ for a continued fraction block\n        base a, b = base(1), c = base(1), d\
    \ = base(0);\n        linfrac(base a): a(a) {}\n        linfrac(base a, base b,\
    \ base c, base d): a(a), b(b), c(c), d(d) {}\n        \n        // composition\
    \ of two linfracs\n        linfrac operator *(linfrac const& t) {\n          \
    \  auto [A, C] = apply(t.a, t.c);\n            auto [B, D] = apply(t.b, t.d);\n\
    \            return {A, B, C, D};\n        }\n        \n        linfrac adj()\
    \ {\n            return {d, -b, -c, a};\n        }\n        \n        // apply\
    \ linfrac to A/B\n        auto apply(base A, base B) {\n            return std::pair{a\
    \ * A + b * B, c * A + d * B};\n        }\n    };\n}\n\n#line 6 \"cp-algo/algebra/modular.hpp\"\
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
    \ x.r % m;\n    }\n}\n\n#line 5 \"cp-algo/algebra/matrix.hpp\"\n#include <valarray>\n\
    #line 9 \"cp-algo/algebra/matrix.hpp\"\n#include <vector>\n#include <array>\n\
    namespace cp_algo::algebra {\n    template<int mod>\n    struct matrix {\n   \
    \     using base = modular<mod>;\n        size_t n, m;\n        std::valarray<std::valarray<base>>\
    \ a;\n        matrix(size_t n, size_t m): n(n), m(m), a(std::valarray<base>(m),\
    \ n) {}\n        matrix(std::valarray<std::valarray<base>> a): n(size(a)), m(n\
    \ ? size(a[0]) : 0), a(a) {}\n\n        auto& operator[] (size_t i) {return a[i];}\n\
    \        auto const& operator[] (size_t i) const {return a[i];}\n        auto&\
    \ row(size_t i) {return a[i];}\n        auto const& row(size_t i) const {return\
    \ a[i];}\n\n        matrix operator -() const {return matrix(-a);}\n        matrix&\
    \ operator *=(base t) {for(auto &it: a) it *= t; return *this;}\n        matrix\
    \ operator *(base t) const {return matrix(*this) *= t;}\n\n        void read()\
    \ {\n            for(size_t i = 0; i < n; i++) {\n                for(size_t j\
    \ = 0; j < m; j++) {\n                    std::cin >> (*this)[i][j];\n       \
    \         }\n            }\n        }\n\n        void print() const {\n      \
    \      for(size_t i = 0; i < n; i++) {\n                for(size_t j = 0; j <\
    \ m; j++) {\n                    std::cout << (*this)[i][j] << \" \\n\"[j + 1\
    \ == m];\n                }\n            }\n        }\n\n        static matrix\
    \ eye(size_t n) {\n            matrix res(n, n);\n            for(size_t i = 0;\
    \ i < n; i++) {\n                res[i][i] = 1;\n            }\n            return\
    \ res;\n        }\n\n        // concatenate matrices\n        matrix operator\
    \ |(matrix const& b) const {\n            assert(n == b.n);\n            matrix\
    \ res(n, m+b.m);\n            for(size_t i = 0; i < n; i++) {\n              \
    \  res[i][std::slice(0,m,1)] = a[i];\n                res[i][std::slice(m,b.m,1)]\
    \ = b[i];\n            }\n            return res;\n        }\n        matrix submatrix(auto\
    \ slicex, auto slicey) const {\n            std::valarray res = a[slicex];\n \
    \           for(auto &row: res) {\n                row = std::valarray(row[slicey]);\n\
    \            }\n            return res;\n        }\n\n        matrix T() const\
    \ {\n            matrix res(m, n);\n            for(size_t i = 0; i < n; i++)\
    \ {\n                for(size_t j = 0; j < m; j++) {\n                    res[j][i]\
    \ = (*this)[i][j];\n                }\n            }\n            return res;\n\
    \        }\n\n        matrix operator *(matrix const& b) const {\n           \
    \ assert(m == b.n);\n            matrix res(n, b.m);\n            for(size_t i\
    \ = 0; i < n; i++) {\n                for(size_t j = 0; j < m; j++) {\n      \
    \              for(size_t k = 0; k < b.m; k++) {\n                        res[i][k].add_unsafe(a[i][j].r\
    \ * b[j][k].r);\n                    }\n                }\n            }\n   \
    \         return res.normalize();\n        }\n\n        matrix pow(uint64_t k)\
    \ const {\n            assert(n == m);\n            return bpow(*this, k, eye(n));\n\
    \        }\n\n        static auto& normalize(auto &a) {\n            for(auto\
    \ &it: a) {\n                it.normalize();\n            }\n            return\
    \ a;\n        }\n        matrix& normalize() {\n            for(auto &it: a) {\n\
    \                normalize(it);\n            }\n            return *this;\n  \
    \      }\n\n        inline static void add_scaled(auto &a, auto const& b, base\
    \ scale, size_t i = 0) {\n            size_t m = size(a);\n            for(; i\
    \ < m; i++) {\n                a[i].add_unsafe(scale.r * b[i].r);\n          \
    \  }\n        }\n\n        enum Mode {normal, reverse};\n        template<Mode\
    \ mode = normal>\n        auto gauss(size_t lim) {\n            size_t rk = 0;\n\
    \            std::vector<size_t> free, pivots;\n            for(size_t i = 0;\
    \ i < lim; i++) {\n                for(size_t j = rk; j < n && a[rk][i].normalize()\
    \ == 0; j++) {\n                    if(a[j][i].normalize() != 0) {\n         \
    \               a[rk] += a[j];\n                    }\n                }\n   \
    \             if(rk == n || normalize(a[rk])[i] == 0) {\n                    free.push_back(i);\n\
    \                } else {\n                    pivots.push_back(i);\n        \
    \            base dinv = -a[rk][i].inv();\n                    for(size_t j =\
    \ mode == reverse ? 0 : rk; j < n; j++) {\n                        if(j != rk)\
    \ {\n                            add_scaled(a[j], a[rk], a[j][i].normalize() *\
    \ dinv, i);\n                        }\n                    }\n              \
    \      rk += 1;\n                }\n            }\n            normalize();\n\
    \            return std::array{pivots, free};\n        }\n        template<Mode\
    \ mode = normal>\n        auto gauss() {\n            return gauss<mode>(m);\n\
    \        }\n\n        size_t rank() const {\n            if(n < m) {\n       \
    \         return T().rank();\n            }\n            return size(matrix(*this).gauss()[0]);\n\
    \        }\n\n        base det() const {\n            assert(n == m);\n      \
    \      matrix b = *this;\n            b.gauss();\n            base res = 1;\n\
    \            for(size_t i = 0; i < n; i++) {\n                res *= b[i][i];\n\
    \            }\n            return res;\n        }\n\n        std::optional<matrix>\
    \ inv() const {\n            assert(n == m);\n            matrix b = *this | eye(n);\n\
    \            if(size(b.gauss<reverse>(n)[0]) < n) {\n                return std::nullopt;\n\
    \            }\n            for(size_t i = 0; i < n; i++) {\n                b[i]\
    \ *= b[i][i].inv();\n            }\n            return b.submatrix(std::slice(0,\
    \ n, 1), std::slice(n, n, 1));\n        }\n\n        // [solution, basis], transposed\n\
    \        std::optional<std::array<matrix, 2>> solve(matrix t) const {\n      \
    \      assert(n == t.n);\n            matrix b = *this | t;\n            auto\
    \ [pivots, free] = b.gauss<reverse>();\n            if(!empty(pivots) && pivots.back()\
    \ >= m) {\n                return std::nullopt;\n            }\n            matrix\
    \ sols(size(free), m);\n            for(size_t j = 0; j < size(pivots); j++) {\n\
    \                base scale = b[j][pivots[j]].inv();\n                for(size_t\
    \ i = 0; i < size(free); i++) {\n                    sols[i][pivots[j]] = b[j][free[i]]\
    \ * scale;\n                }\n            }\n            for(size_t i = 0; free[i]\
    \ < m; i++) {\n                sols[i][free[i]] = -1;\n            }\n       \
    \     return std::array{\n                sols.submatrix(std::slice(size(free)\
    \ - t.m, t.m, 1), std::slice(0, m, 1)),\n                sols.submatrix(std::slice(0,\
    \ size(free) - t.m, 1), std::slice(0, m, 1))\n            };\n        }\n    };\n\
    }\n\n#line 4 \"verify/algebra/matrix/matrix_pow.test.cpp\"\n#include <bits/stdc++.h>\n\
    \nusing namespace std;\nusing namespace cp_algo::algebra;\n\nconst int mod = 998244353;\n\
    \nvoid solve() {\n    int n;\n    uint64_t k;\n    cin >> n >> k;\n    matrix<mod>\
    \ a(n, n);\n    a.read();\n    a.pow(k).print();\n}\n\nsigned main() {\n    //freopen(\"\
    input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n    cin.tie(0);\n \
    \   int t = 1;\n    while(t--) {\n        solve();\n    }\n}\n"
  code: "// @brief Pow of Matrix\n#define PROBLEM \"https://judge.yosupo.jp/problem/pow_of_matrix\"\
    \n#include \"cp-algo/algebra/matrix.hpp\"\n#include <bits/stdc++.h>\n\nusing namespace\
    \ std;\nusing namespace cp_algo::algebra;\n\nconst int mod = 998244353;\n\nvoid\
    \ solve() {\n    int n;\n    uint64_t k;\n    cin >> n >> k;\n    matrix<mod>\
    \ a(n, n);\n    a.read();\n    a.pow(k).print();\n}\n\nsigned main() {\n    //freopen(\"\
    input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n    cin.tie(0);\n \
    \   int t = 1;\n    while(t--) {\n        solve();\n    }\n}"
  dependsOn:
  - cp-algo/algebra/matrix.hpp
  - cp-algo/algebra/common.hpp
  - cp-algo/algebra/modular.hpp
  - cp-algo/random/rng.hpp
  - cp-algo/algebra/affine.hpp
  isVerificationFile: true
  path: verify/algebra/matrix/matrix_pow.test.cpp
  requiredBy: []
  timestamp: '2024-02-11 12:35:24+01:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/algebra/matrix/matrix_pow.test.cpp
layout: document
redirect_from:
- /verify/verify/algebra/matrix/matrix_pow.test.cpp
- /verify/verify/algebra/matrix/matrix_pow.test.cpp.html
title: Pow of Matrix
---
