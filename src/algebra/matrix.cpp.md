---
data:
  _extendedDependsOn: []
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':warning:'
  attributes:
    links: []
  bundledCode: "#line 1 \"src/algebra/matrix.cpp\"\n\nnamespace algebra { // matrix\n\
    \    template<int mod>\n    struct matrix {\n        using base = modular<mod>;\n\
    \        size_t n, m;\n        valarray<valarray<base>> a;\n        matrix(size_t\
    \ n, size_t m): n(n), m(m), a(valarray<base>(m), n) {}\n        matrix(valarray<valarray<base>>\
    \ a): n(size(a)), m(n ? size(a[0]) : 0), a(a) {}\n\n        auto& operator[] (size_t\
    \ i) {return a[i];}\n        auto const& operator[] (size_t i) const {return a[i];}\n\
    \        auto& row(size_t i) {return a[i];}\n        auto const& row(size_t i)\
    \ const {return a[i];}\n\n        matrix operator -() const {return matrix(-a);}\n\
    \        matrix& operator *=(base t) {for(auto &it: a) it *= t; return *this;}\n\
    \        matrix operator *(base t) const {return matrix(*this) *= t;}\n\n    \
    \    void read() {\n            for(size_t i = 0; i < n; i++) {\n            \
    \    for(size_t j = 0; j < m; j++) {\n                    cin >> (*this)[i][j];\n\
    \                }\n            }\n        }\n\n        void print() const {\n\
    \            for(size_t i = 0; i < n; i++) {\n                for(size_t j = 0;\
    \ j < m; j++) {\n                    cout << (*this)[i][j] << \" \\n\"[j + 1 ==\
    \ m];\n                }\n            }\n        }\n\n        static matrix eye(size_t\
    \ n) {\n            matrix res(n, n);\n            for(size_t i = 0; i < n; i++)\
    \ {\n                res[i][i] = 1;\n            }\n            return res;\n\
    \        }\n\n        // concatenate matrices\n        matrix operator |(matrix\
    \ const& b) const {\n            assert(n == b.n);\n            matrix res(n,\
    \ m+b.m);\n            for(size_t i = 0; i < n; i++) {\n                res[i][slice(0,m,1)]\
    \ = a[i];\n                res[i][slice(m,b.m,1)] = b[i];\n            }\n   \
    \         return res;\n        }\n        matrix submatrix(auto slicex, auto slicey)\
    \ const {\n            valarray res = a[slicex];\n            for(auto &row: res)\
    \ {\n                row = valarray(row[slicey]);\n            }\n           \
    \ return res;\n        }\n\n        matrix T() const {\n            matrix res(m,\
    \ n);\n            for(size_t i = 0; i < n; i++) {\n                for(size_t\
    \ j = 0; j < m; j++) {\n                    res[j][i] = (*this)[i][j];\n     \
    \           }\n            }\n            return res;\n        }\n\n        matrix\
    \ operator *(matrix const& b) const {\n            assert(m == b.n);\n       \
    \     matrix res(n, b.m);\n            for(size_t i = 0; i < n; i++) {\n     \
    \           for(size_t j = 0; j < m; j++) {\n                    for(size_t k\
    \ = 0; k < b.m; k++) {\n                        res[i][k].add_unsafe(a[i][j].r\
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
    \            vector<size_t> free, pivots;\n            for(size_t i = 0; i < lim;\
    \ i++) {\n                for(size_t j = rk; j < n && a[rk][i].normalize() ==\
    \ 0; j++) {\n                    if(a[j][i].normalize() != 0) {\n            \
    \            a[rk] += a[j];\n                    }\n                }\n      \
    \          if(rk == n || normalize(a[rk])[i] == 0) {\n                    free.push_back(i);\n\
    \                } else {\n                    pivots.push_back(i);\n        \
    \            base dinv = -a[rk][i].inv();\n                    for(size_t j =\
    \ mode == reverse ? 0 : rk; j < n; j++) {\n                        if(j != rk)\
    \ {\n                            add_scaled(a[j], a[rk], a[j][i].normalize() *\
    \ dinv, i);\n                        }\n                    }\n              \
    \      rk += 1;\n                }\n            }\n            normalize();\n\
    \            return array{pivots, free};\n        }\n        template<Mode mode\
    \ = normal>\n        auto gauss() {\n            return gauss<mode>(m);\n    \
    \    }\n\n        size_t rank() const {\n            if(n < m) {\n           \
    \     return T().rank();\n            }\n            return size(matrix(*this).gauss()[0]);\n\
    \        }\n\n        base det() const {\n            assert(n == m);\n      \
    \      matrix b = *this;\n            b.gauss();\n            base res = 1;\n\
    \            for(size_t i = 0; i < n; i++) {\n                res *= b[i][i];\n\
    \            }\n            return res;\n        }\n\n        optional<matrix>\
    \ inv() const {\n            assert(n == m);\n            matrix b = *this | eye(n);\n\
    \            if(size(b.gauss<reverse>(n)[0]) < n) {\n                return nullopt;\n\
    \            }\n            for(size_t i = 0; i < n; i++) {\n                b[i]\
    \ *= b[i][i].inv();\n            }\n            return b.submatrix(slice(0, n,\
    \ 1), slice(n, n, 1));\n        }\n\n        // [solution, basis], transposed\n\
    \        optional<array<matrix, 2>> solve(matrix t) const {\n            assert(n\
    \ == t.n);\n            matrix b = *this | t;\n            auto [pivots, free]\
    \ = b.gauss<reverse>();\n            if(!empty(pivots) && pivots.back() >= m)\
    \ {\n                return nullopt;\n            }\n            matrix sols(size(free),\
    \ m);\n            for(size_t j = 0; j < size(pivots); j++) {\n              \
    \  base scale = b[j][pivots[j]].inv();\n                for(size_t i = 0; i <\
    \ size(free); i++) {\n                    sols[i][pivots[j]] = b[j][free[i]] *\
    \ scale;\n                }\n            }\n            for(size_t i = 0; free[i]\
    \ < m; i++) {\n                sols[i][free[i]] = -1;\n            }\n       \
    \     return array{\n                sols.submatrix(slice(size(free) - t.m, t.m,\
    \ 1), slice(0, m, 1)),\n                sols.submatrix(slice(0, size(free) - t.m,\
    \ 1), slice(0, m, 1))\n            };\n        }\n    };\n}\n"
  code: "\nnamespace algebra { // matrix\n    template<int mod>\n    struct matrix\
    \ {\n        using base = modular<mod>;\n        size_t n, m;\n        valarray<valarray<base>>\
    \ a;\n        matrix(size_t n, size_t m): n(n), m(m), a(valarray<base>(m), n)\
    \ {}\n        matrix(valarray<valarray<base>> a): n(size(a)), m(n ? size(a[0])\
    \ : 0), a(a) {}\n\n        auto& operator[] (size_t i) {return a[i];}\n      \
    \  auto const& operator[] (size_t i) const {return a[i];}\n        auto& row(size_t\
    \ i) {return a[i];}\n        auto const& row(size_t i) const {return a[i];}\n\n\
    \        matrix operator -() const {return matrix(-a);}\n        matrix& operator\
    \ *=(base t) {for(auto &it: a) it *= t; return *this;}\n        matrix operator\
    \ *(base t) const {return matrix(*this) *= t;}\n\n        void read() {\n    \
    \        for(size_t i = 0; i < n; i++) {\n                for(size_t j = 0; j\
    \ < m; j++) {\n                    cin >> (*this)[i][j];\n                }\n\
    \            }\n        }\n\n        void print() const {\n            for(size_t\
    \ i = 0; i < n; i++) {\n                for(size_t j = 0; j < m; j++) {\n    \
    \                cout << (*this)[i][j] << \" \\n\"[j + 1 == m];\n            \
    \    }\n            }\n        }\n\n        static matrix eye(size_t n) {\n  \
    \          matrix res(n, n);\n            for(size_t i = 0; i < n; i++) {\n  \
    \              res[i][i] = 1;\n            }\n            return res;\n      \
    \  }\n\n        // concatenate matrices\n        matrix operator |(matrix const&\
    \ b) const {\n            assert(n == b.n);\n            matrix res(n, m+b.m);\n\
    \            for(size_t i = 0; i < n; i++) {\n                res[i][slice(0,m,1)]\
    \ = a[i];\n                res[i][slice(m,b.m,1)] = b[i];\n            }\n   \
    \         return res;\n        }\n        matrix submatrix(auto slicex, auto slicey)\
    \ const {\n            valarray res = a[slicex];\n            for(auto &row: res)\
    \ {\n                row = valarray(row[slicey]);\n            }\n           \
    \ return res;\n        }\n\n        matrix T() const {\n            matrix res(m,\
    \ n);\n            for(size_t i = 0; i < n; i++) {\n                for(size_t\
    \ j = 0; j < m; j++) {\n                    res[j][i] = (*this)[i][j];\n     \
    \           }\n            }\n            return res;\n        }\n\n        matrix\
    \ operator *(matrix const& b) const {\n            assert(m == b.n);\n       \
    \     matrix res(n, b.m);\n            for(size_t i = 0; i < n; i++) {\n     \
    \           for(size_t j = 0; j < m; j++) {\n                    for(size_t k\
    \ = 0; k < b.m; k++) {\n                        res[i][k].add_unsafe(a[i][j].r\
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
    \            vector<size_t> free, pivots;\n            for(size_t i = 0; i < lim;\
    \ i++) {\n                for(size_t j = rk; j < n && a[rk][i].normalize() ==\
    \ 0; j++) {\n                    if(a[j][i].normalize() != 0) {\n            \
    \            a[rk] += a[j];\n                    }\n                }\n      \
    \          if(rk == n || normalize(a[rk])[i] == 0) {\n                    free.push_back(i);\n\
    \                } else {\n                    pivots.push_back(i);\n        \
    \            base dinv = -a[rk][i].inv();\n                    for(size_t j =\
    \ mode == reverse ? 0 : rk; j < n; j++) {\n                        if(j != rk)\
    \ {\n                            add_scaled(a[j], a[rk], a[j][i].normalize() *\
    \ dinv, i);\n                        }\n                    }\n              \
    \      rk += 1;\n                }\n            }\n            normalize();\n\
    \            return array{pivots, free};\n        }\n        template<Mode mode\
    \ = normal>\n        auto gauss() {\n            return gauss<mode>(m);\n    \
    \    }\n\n        size_t rank() const {\n            if(n < m) {\n           \
    \     return T().rank();\n            }\n            return size(matrix(*this).gauss()[0]);\n\
    \        }\n\n        base det() const {\n            assert(n == m);\n      \
    \      matrix b = *this;\n            b.gauss();\n            base res = 1;\n\
    \            for(size_t i = 0; i < n; i++) {\n                res *= b[i][i];\n\
    \            }\n            return res;\n        }\n\n        optional<matrix>\
    \ inv() const {\n            assert(n == m);\n            matrix b = *this | eye(n);\n\
    \            if(size(b.gauss<reverse>(n)[0]) < n) {\n                return nullopt;\n\
    \            }\n            for(size_t i = 0; i < n; i++) {\n                b[i]\
    \ *= b[i][i].inv();\n            }\n            return b.submatrix(slice(0, n,\
    \ 1), slice(n, n, 1));\n        }\n\n        // [solution, basis], transposed\n\
    \        optional<array<matrix, 2>> solve(matrix t) const {\n            assert(n\
    \ == t.n);\n            matrix b = *this | t;\n            auto [pivots, free]\
    \ = b.gauss<reverse>();\n            if(!empty(pivots) && pivots.back() >= m)\
    \ {\n                return nullopt;\n            }\n            matrix sols(size(free),\
    \ m);\n            for(size_t j = 0; j < size(pivots); j++) {\n              \
    \  base scale = b[j][pivots[j]].inv();\n                for(size_t i = 0; i <\
    \ size(free); i++) {\n                    sols[i][pivots[j]] = b[j][free[i]] *\
    \ scale;\n                }\n            }\n            for(size_t i = 0; free[i]\
    \ < m; i++) {\n                sols[i][free[i]] = -1;\n            }\n       \
    \     return array{\n                sols.submatrix(slice(size(free) - t.m, t.m,\
    \ 1), slice(0, m, 1)),\n                sols.submatrix(slice(0, size(free) - t.m,\
    \ 1), slice(0, m, 1))\n            };\n        }\n    };\n}\n"
  dependsOn: []
  isVerificationFile: false
  path: src/algebra/matrix.cpp
  requiredBy: []
  timestamp: '2024-02-10 16:42:18+01:00'
  verificationStatus: LIBRARY_NO_TESTS
  verifiedWith: []
documentation_of: src/algebra/matrix.cpp
layout: document
redirect_from:
- /library/src/algebra/matrix.cpp
- /library/src/algebra/matrix.cpp.html
title: src/algebra/matrix.cpp
---
