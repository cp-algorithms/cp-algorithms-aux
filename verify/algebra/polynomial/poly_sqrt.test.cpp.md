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
    path: cp-algo/algebra/fft.hpp
    title: cp-algo/algebra/fft.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/polynomial.hpp
    title: cp-algo/algebra/polynomial.hpp
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
    PROBLEM: https://judge.yosupo.jp/problem/sqrt_of_formal_power_series
    document_title: Sqrt of Formal Power Series
    links:
    - https://judge.yosupo.jp/problem/sqrt_of_formal_power_series
  bundledCode: "#line 1 \"verify/algebra/polynomial/poly_sqrt.test.cpp\"\n// @brief\
    \ Sqrt of Formal Power Series\n#define PROBLEM \"https://judge.yosupo.jp/problem/sqrt_of_formal_power_series\"\
    \n#line 1 \"cp-algo/algebra/polynomial.hpp\"\n\n\n#line 1 \"cp-algo/algebra/common.hpp\"\
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
    \     init = true;\n        }\n        return F[n];\n    }\n}\n\n#line 1 \"cp-algo/algebra/fft.hpp\"\
    \n\n\n#line 1 \"cp-algo/algebra/modular.hpp\"\n\n\n#line 1 \"cp-algo/random/rng.hpp\"\
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
    \ {\n            return a * x + b;\n        }\n    };\n}\n\n#line 6 \"cp-algo/algebra/modular.hpp\"\
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
    \ x.r % m;\n    }\n}\n\n#line 7 \"cp-algo/algebra/fft.hpp\"\n#include <vector>\n\
    namespace cp_algo::algebra::fft {\n    using ftype = double;\n    struct point\
    \ {\n        ftype x, y;\n        \n        ftype real() {return x;}\n       \
    \ ftype imag() {return y;}\n        \n        point(): x(0), y(0){}\n        point(ftype\
    \ x, ftype y = 0): x(x), y(y){}\n        \n        static point polar(ftype rho,\
    \ ftype ang) {\n            return point{rho * cos(ang), rho * sin(ang)};\n  \
    \      }\n        \n        point conj() const {\n            return {x, -y};\n\
    \        }\n        \n        point operator +=(const point &t) {x += t.x, y +=\
    \ t.y; return *this;}\n        point operator +(const point &t) const {return\
    \ point(*this) += t;}\n        point operator -(const point &t) const {return\
    \ {x - t.x, y - t.y};}\n        point operator *(const point &t) const {return\
    \ {x * t.x - y * t.y, x * t.y + y * t.x};}\n    };\n\n    point w[maxn]; // w[2^n\
    \ + k] = exp(pi * k / (2^n))\n    int bitr[maxn];// b[2^n + k] = bitreverse(k)\n\
    \    const ftype pi = acos(-1);\n    bool initiated = 0;\n    void init() {\n\
    \        if(!initiated) {\n            for(int i = 1; i < maxn; i *= 2) {\n  \
    \              int ti = i / 2;\n                for(int j = 0; j < i; j++) {\n\
    \                    w[i + j] = point::polar(ftype(1), pi * j / i);\n        \
    \            if(ti) {\n                        bitr[i + j] = 2 * bitr[ti + j %\
    \ ti] + (j >= ti);\n                    }\n                }\n            }\n\
    \            initiated = 1;\n        }\n    }\n    \n    void fft(auto &a, int\
    \ n) {\n        init();\n        if(n == 1) {\n            return;\n        }\n\
    \        int hn = n / 2;\n        for(int i = 0; i < n; i++) {\n            int\
    \ ti = 2 * bitr[hn + i % hn] + (i > hn);\n            if(i < ti) {\n         \
    \       std::swap(a[i], a[ti]);\n            }\n        }\n        for(int i =\
    \ 1; i < n; i *= 2) {\n            for(int j = 0; j < n; j += 2 * i) {\n     \
    \           for(int k = j; k < j + i; k++) {\n                    point t = a[k\
    \ + i] * w[i + k - j];\n                    a[k + i] = a[k] - t;\n           \
    \         a[k] += t;\n                }\n            }\n        }\n    }\n   \
    \ \n    void mul_slow(std::vector<auto> &a, const std::vector<auto> &b) {\n  \
    \      if(a.empty() || b.empty()) {\n            a.clear();\n        } else {\n\
    \            int n = a.size();\n            int m = b.size();\n            a.resize(n\
    \ + m - 1);\n            for(int k = n + m - 2; k >= 0; k--) {\n             \
    \   a[k] *= b[0];\n                for(int j = std::max(k - n + 1, 1); j < std::min(k\
    \ + 1, m); j++) {\n                    a[k] += a[k - j] * b[j];\n            \
    \    }\n            }\n        }\n    }\n    \n    template<int m>\n    struct\
    \ dft {\n        static constexpr int split = 1 << 15;\n        std::vector<point>\
    \ A;\n        \n        dft(std::vector<modular<m>> const& a, size_t n): A(n)\
    \ {\n            for(size_t i = 0; i < std::min(n, a.size()); i++) {\n       \
    \         A[i] = point(\n                    a[i].rem() % split,\n           \
    \         a[i].rem() / split\n                );\n            }\n            if(n)\
    \ {\n                fft(A, n);\n            }\n        }\n    \n        auto\
    \ operator * (dft const& B) {\n            assert(A.size() == B.A.size());\n \
    \           size_t n = A.size();\n            if(!n) {\n                return\
    \ std::vector<modular<m>>();\n            }\n            std::vector<point> C(n),\
    \ D(n);\n            for(size_t i = 0; i < n; i++) {\n                C[i] = A[i]\
    \ * (B[i] + B[(n - i) % n].conj());\n                D[i] = A[i] * (B[i] - B[(n\
    \ - i) % n].conj());\n            }\n            fft(C, n);\n            fft(D,\
    \ n);\n            reverse(begin(C) + 1, end(C));\n            reverse(begin(D)\
    \ + 1, end(D));\n            int t = 2 * n;\n            std::vector<modular<m>>\
    \ res(n);\n            for(size_t i = 0; i < n; i++) {\n                modular<m>\
    \ A0 = llround(C[i].real() / t);\n                modular<m> A1 = llround(C[i].imag()\
    \ / t + D[i].imag() / t);\n                modular<m> A2 = llround(D[i].real()\
    \ / t);\n                res[i] = A0 + A1 * split - A2 * split * split;\n    \
    \        }\n            return res;\n        }\n        \n        point& operator\
    \ [](int i) {return A[i];}\n        point operator [](int i) const {return A[i];}\n\
    \    };\n    \n    size_t com_size(size_t as, size_t bs) {\n        if(!as ||\
    \ !bs) {\n            return 0;\n        }\n        size_t n = as + bs - 1;\n\
    \        while(__builtin_popcount(n) != 1) {\n            n++;\n        }\n  \
    \      return n;\n    }\n    \n    template<int m>\n    void mul(std::vector<modular<m>>\
    \ &a, std::vector<modular<m>> b) {\n        if(std::min(a.size(), b.size()) <\
    \ magic) {\n            mul_slow(a, b);\n            return;\n        }\n    \
    \    auto n = com_size(a.size(), b.size());\n        auto A = dft<m>(a, n);\n\
    \        if(a == b) {\n            a = A * A;\n        } else {\n            a\
    \ = A * dft<m>(b, n);\n        }\n    }\n}\n\n#line 5 \"cp-algo/algebra/polynomial.hpp\"\
    \n#include <functional>\n#line 9 \"cp-algo/algebra/polynomial.hpp\"\n#include\
    \ <utility>\n#line 11 \"cp-algo/algebra/polynomial.hpp\"\nnamespace cp_algo::algebra\
    \ {\n    template<typename T>\n    struct poly {\n        std::vector<T> a;\n\
    \        \n        void normalize() { // get rid of leading zeroes\n         \
    \   while(!a.empty() && a.back() == T(0)) {\n                a.pop_back();\n \
    \           }\n        }\n        \n        poly(){}\n        poly(T a0) : a{a0}{normalize();}\n\
    \        poly(const std::vector<T> &t) : a(t){normalize();}\n        \n      \
    \  poly operator -() const {\n            auto t = *this;\n            for(auto\
    \ &it: t.a) {\n                it = -it;\n            }\n            return t;\n\
    \        }\n        \n        poly operator += (const poly &t) {\n           \
    \ a.resize(std::max(a.size(), t.a.size()));\n            for(size_t i = 0; i <\
    \ t.a.size(); i++) {\n                a[i] += t.a[i];\n            }\n       \
    \     normalize();\n            return *this;\n        }\n        \n        poly\
    \ operator -= (const poly &t) {\n            a.resize(std::max(a.size(), t.a.size()));\n\
    \            for(size_t i = 0; i < t.a.size(); i++) {\n                a[i] -=\
    \ t.a[i];\n            }\n            normalize();\n            return *this;\n\
    \        }\n        poly operator + (const poly &t) const {return poly(*this)\
    \ += t;}\n        poly operator - (const poly &t) const {return poly(*this) -=\
    \ t;}\n        \n        poly mod_xk(size_t k) const { // get first k coefficients\n\
    \            return std::vector<T>(begin(a), begin(a) + std::min(k, a.size()));\n\
    \        }\n        \n        poly mul_xk(size_t k) const { // multiply by x^k\n\
    \            auto res = a;\n            res.insert(begin(res), k, 0);\n      \
    \      return res;\n        }\n        \n        poly div_xk(size_t k) const {\
    \ // drop first k coefficients\n            return std::vector<T>(begin(a) + std::min(k,\
    \ a.size()), end(a));\n        }\n        \n        poly substr(size_t l, size_t\
    \ r) const { // return mod_xk(r).div_xk(l)\n            return std::vector<T>(\n\
    \                begin(a) + std::min(l, a.size()),\n                begin(a) +\
    \ std::min(r, a.size())\n            );\n        }\n        \n        poly operator\
    \ *= (const poly &t) {fft::mul(a, t.a); normalize(); return *this;}\n        poly\
    \ operator * (const poly &t) const {return poly(*this) *= t;}\n        \n    \
    \    poly reverse(size_t n) const { // computes x^n A(x^{-1})\n            auto\
    \ res = a;\n            res.resize(std::max(n, res.size()));\n            return\
    \ std::vector<T>(res.rbegin(), res.rbegin() + n);\n        }\n        \n     \
    \   poly reverse() const {\n            return reverse(deg() + 1);\n        }\n\
    \        \n        std::pair<poly, poly> divmod_slow(const poly &b) const { //\
    \ when divisor or quotient is small\n            std::vector<T> A(a);\n      \
    \      std::vector<T> res;\n            T b_lead_inv = b.a.back().inv();\n   \
    \         while(A.size() >= b.a.size()) {\n                res.push_back(A.back()\
    \ * b_lead_inv);\n                if(res.back() != T(0)) {\n                 \
    \   for(size_t i = 0; i < b.a.size(); i++) {\n                        A[A.size()\
    \ - i - 1] -= res.back() * b.a[b.a.size() - i - 1];\n                    }\n \
    \               }\n                A.pop_back();\n            }\n            std::reverse(begin(res),\
    \ end(res));\n            return {res, A};\n        }\n        \n        std::pair<poly,\
    \ poly> divmod_hint(poly const& b, poly const& binv) const { // when inverse is\
    \ known\n            assert(!b.is_zero());\n            if(deg() < b.deg()) {\n\
    \                return {poly{0}, *this};\n            }\n            int d =\
    \ deg() - b.deg();\n            if(std::min(d, b.deg()) < magic) {\n         \
    \       return divmod_slow(b);\n            }\n            poly D = (reverse().mod_xk(d\
    \ + 1) * binv.mod_xk(d + 1)).mod_xk(d + 1).reverse(d + 1);\n            return\
    \ {D, *this - D * b};\n        }\n        \n        std::pair<poly, poly> divmod(const\
    \ poly &b) const { // returns quotiend and remainder of a mod b\n            assert(!b.is_zero());\n\
    \            if(deg() < b.deg()) {\n                return {poly{0}, *this};\n\
    \            }\n            int d = deg() - b.deg();\n            if(std::min(d,\
    \ b.deg()) < magic) {\n                return divmod_slow(b);\n            }\n\
    \            poly D = (reverse().mod_xk(d + 1) * b.reverse().inv(d + 1)).mod_xk(d\
    \ + 1).reverse(d + 1);\n            return {D, *this - D * b};\n        }\n  \
    \      \n        // (ax+b) / (cx+d)\n        struct transform {\n            poly\
    \ a, b, c, d;\n            transform(poly a, poly b = T(1), poly c = T(1), poly\
    \ d = T(0)): a(a), b(b), c(c), d(d){}\n            \n            transform operator\
    \ *(transform const& t) {\n                return {\n                    a*t.a\
    \ + b*t.c, a*t.b + b*t.d,\n                    c*t.a + d*t.c, c*t.b + d*t.d\n\
    \                };\n            }\n            \n            transform adj()\
    \ {\n                return transform(d, -b, -c, a);\n            }\n        \
    \    \n            auto apply(poly A, poly B) {\n                return make_pair(a\
    \ * A + b * B, c * A + d * B);\n            }\n        };\n        \n        template<typename\
    \ Q>\n        static void concat(std::vector<Q> &a, std::vector<Q> const& b) {\n\
    \            for(auto it: b) {\n                a.push_back(it);\n           \
    \ }\n        }\n        \n        // finds a transform that changes A/B to A'/B'\
    \ such that\n        // deg B' is at least 2 times less than deg A\n        static\
    \ std::pair<std::vector<poly>, transform> half_gcd(poly A, poly B) {\n       \
    \     assert(A.deg() >= B.deg());\n            int m = (A.deg() + 1) / 2;\n  \
    \          if(B.deg() < m) {\n                return {{}, {T(1), T(0), T(0), T(1)}};\n\
    \            }\n            auto [ar, Tr] = half_gcd(A.div_xk(m), B.div_xk(m));\n\
    \            tie(A, B) = Tr.adj().apply(A, B);\n            if(B.deg() < m) {\n\
    \                return {ar, Tr};\n            }\n            auto [ai, R] = A.divmod(B);\n\
    \            tie(A, B) = make_pair(B, R);\n            int k = 2 * m - B.deg();\n\
    \            auto [as, Ts] = half_gcd(A.div_xk(k), B.div_xk(k));\n           \
    \ concat(ar, {ai});\n            concat(ar, as);\n            return {ar, Tr *\
    \ transform(ai) * Ts};\n        }\n        \n        // return a transform that\
    \ reduces A / B to gcd(A, B) / 0\n        static std::pair<std::vector<poly>,\
    \ transform> full_gcd(poly A, poly B) {\n            std::vector<poly> ak;\n \
    \           std::vector<transform> trs;\n            while(!B.is_zero()) {\n \
    \               if(2 * B.deg() > A.deg()) {\n                    auto [a, Tr]\
    \ = half_gcd(A, B);\n                    concat(ak, a);\n                    trs.push_back(Tr);\n\
    \                    tie(A, B) = trs.back().adj().apply(A, B);\n             \
    \   } else {\n                    auto [a, R] = A.divmod(B);\n               \
    \     ak.push_back(a);\n                    trs.emplace_back(a);\n           \
    \         tie(A, B) = make_pair(B, R);\n                }\n            }\n   \
    \         trs.emplace_back(T(1), T(0), T(0), T(1));\n            while(trs.size()\
    \ >= 2) {\n                trs[trs.size() - 2] = trs[trs.size() - 2] * trs[trs.size()\
    \ - 1];\n                trs.pop_back();\n            }\n            return {ak,\
    \ trs.back()};\n        }\n                \n        static poly gcd(poly A, poly\
    \ B) {\n            if(A.deg() < B.deg()) {\n                return gcd(B, A);\n\
    \            }\n            auto [a, Tr] = full_gcd(A, B);\n            return\
    \ Tr.d * A - Tr.b * B;\n        }\n\n        \n        // Returns the characteristic\
    \ polynomial\n        // of the minimum linear recurrence for the sequence\n \
    \       poly min_rec_slow(int d) const {\n            auto R1 = mod_xk(d + 1).reverse(d\
    \ + 1), R2 = xk(d + 1);\n            auto Q1 = poly(T(1)), Q2 = poly(T(0));\n\
    \            while(!R2.is_zero()) {\n                auto [a, nR] = R1.divmod(R2);\
    \ // R1 = a*R2 + nR, deg nR < deg R2\n                tie(R1, R2) = make_tuple(R2,\
    \ nR);\n                tie(Q1, Q2) = make_tuple(Q2, Q1 + a * Q2);\n         \
    \       if(R2.deg() < Q2.deg()) {\n                    return Q2 / Q2.lead();\n\
    \                }\n            }\n            assert(0);\n        }\n       \
    \ \n        static transform convergent(auto L, auto R) { // computes product\
    \ on [L, R)\n            if(R - L == 1) {\n                return transform(*L);\n\
    \            } else {\n                int s = 0;\n                for(int i =\
    \ 0; i < R - L; i++) {\n                    s += L[i].a.size();\n            \
    \    }\n                int c = 0;\n                for(int i = 0; i < R - L;\
    \ i++) {\n                    c += L[i].a.size();\n                    if(2 *\
    \ c > s) {\n                        return convergent(L, L + i) * convergent(L\
    \ + i, R);\n                    }\n                }\n                assert(0);\n\
    \            }\n        }\n        \n        poly min_rec(int d) const {\n   \
    \         if(d < magic) {\n                return min_rec_slow(d);\n         \
    \   }\n            auto R2 = mod_xk(d + 1).reverse(d + 1), R1 = xk(d + 1);\n \
    \           if(R2.is_zero()) {\n                return poly(1);\n            }\n\
    \            auto [a, Tr] = full_gcd(R1, R2);\n            int dr = (d + 1) -\
    \ a[0].deg();\n            int dp = 0;\n            for(size_t i = 0; i + 1 <\
    \ a.size(); i++) {\n                dr -= a[i + 1].deg();\n                dp\
    \ += a[i].deg();\n                if(dr < dp) {\n                    auto ans\
    \ = convergent(begin(a), begin(a) + i + 1);\n                    return ans.a\
    \ / ans.a.lead();\n                }\n            }\n            auto ans = convergent(begin(a),\
    \ end(a));\n            return ans.a / ans.a.lead();\n        }\n        \n  \
    \      // calculate inv to *this modulo t\n        // quadratic complexity\n \
    \       std::optional<poly> inv_mod_slow(poly const& t) const {\n            auto\
    \ R1 = *this, R2 = t;\n            auto Q1 = poly(T(1)), Q2 = poly(T(0));\n  \
    \          int k = 0;\n            while(!R2.is_zero()) {\n                k ^=\
    \ 1;\n                auto [a, nR] = R1.divmod(R2);\n                tie(R1, R2)\
    \ = make_tuple(R2, nR);\n                tie(Q1, Q2) = make_tuple(Q2, Q1 + a *\
    \ Q2);\n            }\n            if(R1.deg() > 0) {\n                return\
    \ std::nullopt;\n            } else {\n                return (k ? -Q1 : Q1) /\
    \ R1[0];\n            }\n        }\n        \n        std::optional<poly> inv_mod(poly\
    \ const &t) const {\n            assert(!t.is_zero());\n            if(false &&\
    \ std::min(deg(), t.deg()) < magic) {\n                return inv_mod_slow(t);\n\
    \            }\n            auto A = t, B = *this % t;\n            auto [a, Tr]\
    \ = full_gcd(A, B);\n            auto g = Tr.d * A - Tr.b * B;\n            if(g.deg()\
    \ != 0) {\n                return std::nullopt;\n            }\n            return\
    \ -Tr.b / g[0];\n        };\n        \n        poly operator / (const poly &t)\
    \ const {return divmod(t).first;}\n        poly operator % (const poly &t) const\
    \ {return divmod(t).second;}\n        poly operator /= (const poly &t) {return\
    \ *this = divmod(t).first;}\n        poly operator %= (const poly &t) {return\
    \ *this = divmod(t).second;}\n        poly operator *= (const T &x) {\n      \
    \      for(auto &it: a) {\n                it *= x;\n            }\n         \
    \   normalize();\n            return *this;\n        }\n        poly operator\
    \ /= (const T &x) {\n            return *this *= x.inv();\n        }\n       \
    \ poly operator * (const T &x) const {return poly(*this) *= x;}\n        poly\
    \ operator / (const T &x) const {return poly(*this) /= x;}\n        \n       \
    \ poly conj() const { // A(x) -> A(-x)\n            auto res = *this;\n      \
    \      for(int i = 1; i <= deg(); i += 2) {\n                res.a[i] = -res[i];\n\
    \            }\n            return res;\n        }\n        \n        void print(int\
    \ n) const {\n            for(int i = 0; i < n; i++) {\n                std::cout\
    \ << (*this)[i] << ' ';\n            }\n            std::cout << \"\\n\";\n  \
    \      }\n        \n        void print() const {\n            print(deg() + 1);\n\
    \        }\n        \n        T eval(T x) const { // evaluates in single point\
    \ x\n            T res(0);\n            for(int i = deg(); i >= 0; i--) {\n  \
    \              res *= x;\n                res += a[i];\n            }\n      \
    \      return res;\n        }\n        \n        T lead() const { // leading coefficient\n\
    \            assert(!is_zero());\n            return a.back();\n        }\n  \
    \      \n        int deg() const { // degree, -1 for P(x) = 0\n            return\
    \ (int)a.size() - 1;\n        }\n        \n        bool is_zero() const {\n  \
    \          return a.empty();\n        }\n        \n        T operator [](int idx)\
    \ const {\n            return idx < 0 || idx > deg() ? T(0) : a[idx];\n      \
    \  }\n        \n        T& coef(size_t idx) { // mutable reference at coefficient\n\
    \            return a[idx];\n        }\n        \n        bool operator == (const\
    \ poly &t) const {return a == t.a;}\n        bool operator != (const poly &t)\
    \ const {return a != t.a;}\n        \n        poly deriv(int k = 1) { // calculate\
    \ derivative\n            if(deg() + 1 < k) {\n                return poly(T(0));\n\
    \            }\n            std::vector<T> res(deg() + 1 - k);\n            for(int\
    \ i = k; i <= deg(); i++) {\n                res[i - k] = fact<T>(i) * rfact<T>(i\
    \ - k) * a[i];\n            }\n            return res;\n        }\n        \n\
    \        poly integr() { // calculate integral with C = 0\n            std::vector<T>\
    \ res(deg() + 2);\n            for(int i = 0; i <= deg(); i++) {\n           \
    \     res[i + 1] = a[i] * small_inv<T>(i + 1);\n            }\n            return\
    \ res;\n        }\n        \n        size_t trailing_xk() const { // Let p(x)\
    \ = x^k * t(x), return k\n            if(is_zero()) {\n                return\
    \ -1;\n            }\n            int res = 0;\n            while(a[res] == T(0))\
    \ {\n                res++;\n            }\n            return res;\n        }\n\
    \        \n        poly log(size_t n) { // calculate log p(x) mod x^n\n      \
    \      assert(a[0] == T(1));\n            return (deriv().mod_xk(n) * inv(n)).integr().mod_xk(n);\n\
    \        }\n        \n        poly exp(size_t n) { // calculate exp p(x) mod x^n\n\
    \            if(is_zero()) {\n                return T(1);\n            }\n  \
    \          assert(a[0] == T(0));\n            poly ans = T(1);\n            size_t\
    \ a = 1;\n            while(a < n) {\n                poly C = ans.log(2 * a).div_xk(a)\
    \ - substr(a, 2 * a);\n                ans -= (ans * C).mod_xk(a).mul_xk(a);\n\
    \                a *= 2;\n            }\n            return ans.mod_xk(n);\n \
    \       }\n        \n        poly pow_bin(int64_t k, size_t n) { // O(n log n\
    \ log k)\n            if(k == 0) {\n                return poly(1).mod_xk(n);\n\
    \            } else {\n                auto t = pow(k / 2, n);\n             \
    \   t = (t * t).mod_xk(n);\n                return (k % 2 ? *this * t : t).mod_xk(n);\n\
    \            }\n        }\n        \n        // Do not compute inverse from scratch\n\
    \        poly powmod_hint(int64_t k, poly const& md, poly const& mdinv) {\n  \
    \          if(k == 0) {\n                return poly(1);\n            } else {\n\
    \                auto t = powmod_hint(k / 2, md, mdinv);\n                t =\
    \ (t * t).divmod_hint(md, mdinv).second;\n                if(k % 2) {\n      \
    \              t = (t * *this).divmod_hint(md, mdinv).second;\n              \
    \  }\n                return t;\n            }\n        }\n\n        poly circular_closure(size_t\
    \ m) const {\n            if(deg() == -1) {\n                return *this;\n \
    \           }\n            auto t = *this;\n            for(size_t i = t.deg();\
    \ i >= m; i--) {\n                t.a[i - m] += t.a[i];\n            }\n     \
    \       t.a.resize(std::min(t.a.size(), m));\n            return t;\n        }\n\
    \n        static poly mul_circular(poly const& a, poly const& b, size_t m) {\n\
    \            return (a.circular_closure(m) * b.circular_closure(m)).circular_closure(m);\n\
    \        }\n\n        poly powmod_circular(int64_t k, size_t m) {\n          \
    \  if(k == 0) {\n                return poly(1);\n            } else {\n     \
    \           auto t = powmod_circular(k / 2, m);\n                t = mul_circular(t,\
    \ t, m);\n                if(k % 2) {\n                    t = mul_circular(t,\
    \ *this, m);\n                }\n                return t;\n            }\n  \
    \      }\n        \n        poly powmod(int64_t k, poly const& md) {\n       \
    \     int d = md.deg();\n            if(d == -1) {\n                return k ?\
    \ *this : poly(T(1));\n            }\n            if(md == xk(d)) {\n        \
    \        return pow(k, d);\n            }\n            if(md == xk(d) - poly(T(1)))\
    \ {\n                return powmod_circular(k, d);\n            }\n          \
    \  auto mdinv = md.reverse().inv(md.deg() + 1);\n            return powmod_hint(k,\
    \ md, mdinv);\n        }\n        \n        // O(d * n) with the derivative trick\
    \ from\n        // https://codeforces.com/blog/entry/73947?#comment-581173\n \
    \       poly pow_dn(int64_t k, size_t n) {\n            if(n == 0) {\n       \
    \         return poly(T(0));\n            }\n            assert((*this)[0] !=\
    \ T(0));\n            std::vector<T> Q(n);\n            Q[0] = bpow(a[0], k);\n\
    \            auto a0inv = a[0].inv();\n            for(int i = 1; i < (int)n;\
    \ i++) {\n                for(int j = 1; j <= std::min(deg(), i); j++) {\n   \
    \                 Q[i] += a[j] * Q[i - j] * (T(k) * T(j) - T(i - j));\n      \
    \          }\n                Q[i] *= small_inv<T>(i) * a0inv;\n            }\n\
    \            return Q;\n        }\n        \n        // calculate p^k(n) mod x^n\
    \ in O(n log n)\n        // might be quite slow due to high constant\n       \
    \ poly pow(int64_t k, size_t n) {\n            if(is_zero()) {\n             \
    \   return k ? *this : poly(1);\n            }\n            int i = trailing_xk();\n\
    \            if(i > 0) {\n                return k >= int64_t(n + i - 1) / i ?\
    \ poly(T(0)) : div_xk(i).pow(k, n - i * k).mul_xk(i * k);\n            }\n   \
    \         if(std::min(deg(), (int)n) <= magic) {\n                return pow_dn(k,\
    \ n);\n            }\n            if(k <= magic) {\n                return pow_bin(k,\
    \ n);\n            }\n            T j = a[i];\n            poly t = *this / j;\n\
    \            return bpow(j, k) * (t.log(n) * T(k)).exp(n).mod_xk(n);\n       \
    \ }\n        \n        // returns std::nullopt if undefined\n        std::optional<poly>\
    \ sqrt(size_t n) const {\n            if(is_zero()) {\n                return\
    \ *this;\n            }\n            int i = trailing_xk();\n            if(i\
    \ % 2) {\n                return std::nullopt;\n            } else if(i > 0) {\n\
    \                auto ans = div_xk(i).sqrt(n - i / 2);\n                return\
    \ ans ? ans->mul_xk(i / 2) : ans;\n            }\n            auto st = (*this)[0].sqrt();\n\
    \            if(st) {\n                poly ans = *st;\n                size_t\
    \ a = 1;\n                while(a < n) {\n                    a *= 2;\n      \
    \              ans -= (ans - mod_xk(a) * ans.inv(a)).mod_xk(a) / 2;\n        \
    \        }\n                return ans.mod_xk(n);\n            }\n           \
    \ return std::nullopt;\n        }\n        \n        poly mulx(T a) const { //\
    \ component-wise multiplication with a^k\n            T cur = 1;\n           \
    \ poly res(*this);\n            for(int i = 0; i <= deg(); i++) {\n          \
    \      res.coef(i) *= cur;\n                cur *= a;\n            }\n       \
    \     return res;\n        }\n\n        poly mulx_sq(T a) const { // component-wise\
    \ multiplication with a^{k choose 2}\n            T cur = 1, total = 1;\n    \
    \        poly res(*this);\n            for(int i = 0; i <= deg(); i++) {\n   \
    \             res.coef(i) *= total;\n                cur *= a;\n             \
    \   total *= cur;\n            }\n            return res;\n        }\n\n     \
    \   // be mindful of maxn, as the function\n        // requires multiplying polynomials\
    \ of size deg() and n+deg()!\n        poly chirpz(T z, int n) const { // P(1),\
    \ P(z), P(z^2), ..., P(z^(n-1))\n            if(is_zero()) {\n               \
    \ return std::vector<T>(n, 0);\n            }\n            if(z == T(0)) {\n \
    \               std::vector<T> ans(n, (*this)[0]);\n                if(n > 0)\
    \ {\n                    ans[0] = accumulate(begin(a), end(a), T(0));\n      \
    \          }\n                return ans;\n            }\n            auto A =\
    \ mulx_sq(z.inv());\n            auto B = ones(n+deg()).mulx_sq(z);\n        \
    \    return semicorr(B, A).mod_xk(n).mulx_sq(z.inv());\n        }\n\n        //\
    \ res[i] = prod_{1 <= j <= i} 1/(1 - z^j)\n        static auto _1mzk_prod_inv(T\
    \ z, int n) {\n            std::vector<T> res(n, 1), zk(n);\n            zk[0]\
    \ = 1;\n            for(int i = 1; i < n; i++) {\n                zk[i] = zk[i\
    \ - 1] * z;\n                res[i] = res[i - 1] * (T(1) - zk[i]);\n         \
    \   }\n            res.back() = res.back().inv();\n            for(int i = n -\
    \ 2; i >= 0; i--) {\n                res[i] = (T(1) - zk[i+1]) * res[i+1];\n \
    \           }\n            return res;\n        }\n        \n        // prod_{0\
    \ <= j < n} (1 - z^j x)\n        static auto _1mzkx_prod(T z, int n) {\n     \
    \       if(n == 1) {\n                return poly(std::vector<T>{1, -1});\n  \
    \          } else {\n                auto t = _1mzkx_prod(z, n / 2);\n       \
    \         t *= t.mulx(bpow(z, n / 2));\n                if(n % 2) {\n        \
    \            t *= poly(std::vector<T>{1, -bpow(z, n - 1)});\n                }\n\
    \                return t;\n            }\n        }\n\n        poly chirpz_inverse(T\
    \ z, int n) const { // P(1), P(z), P(z^2), ..., P(z^(n-1))\n            if(is_zero())\
    \ {\n                return {};\n            }\n            if(z == T(0)) {\n\
    \                if(n == 1) {\n                    return *this;\n           \
    \     } else {\n                    return std::vector{(*this)[1], (*this)[0]\
    \ - (*this)[1]};\n                }\n            }\n            std::vector<T>\
    \ y(n);\n            for(int i = 0; i < n; i++) {\n                y[i] = (*this)[i];\n\
    \            }\n            auto prods_pos = _1mzk_prod_inv(z, n);\n         \
    \   auto prods_neg = _1mzk_prod_inv(z.inv(), n);\n\n            T zn = bpow(z,\
    \ n-1).inv();\n            T znk = 1;\n            for(int i = 0; i < n; i++)\
    \ {\n                y[i] *= znk * prods_neg[i] * prods_pos[(n - 1) - i];\n  \
    \              znk *= zn;\n            }\n\n            poly p_over_q = poly(y).chirpz(z,\
    \ n);\n            poly q = _1mzkx_prod(z, n);\n\n            return (p_over_q\
    \ * q).mod_xk(n).reverse(n);\n        }\n\n        static poly build(std::vector<poly>\
    \ &res, int v, auto L, auto R) { // builds evaluation tree for (x-a1)(x-a2)...(x-an)\n\
    \            if(R - L == 1) {\n                return res[v] = std::vector<T>{-*L,\
    \ 1};\n            } else {\n                auto M = L + (R - L) / 2;\n     \
    \           return res[v] = build(res, 2 * v, L, M) * build(res, 2 * v + 1, M,\
    \ R);\n            }\n        }\n\n        poly to_newton(std::vector<poly> &tree,\
    \ int v, auto l, auto r) {\n            if(r - l == 1) {\n                return\
    \ *this;\n            } else {\n                auto m = l + (r - l) / 2;\n  \
    \              auto A = (*this % tree[2 * v]).to_newton(tree, 2 * v, l, m);\n\
    \                auto B = (*this / tree[2 * v]).to_newton(tree, 2 * v + 1, m,\
    \ r);\n                return A + B.mul_xk(m - l);\n            }\n        }\n\
    \n        poly to_newton(std::vector<T> p) {\n            if(is_zero()) {\n  \
    \              return *this;\n            }\n            int n = p.size();\n \
    \           std::vector<poly> tree(4 * n);\n            build(tree, 1, begin(p),\
    \ end(p));\n            return to_newton(tree, 1, begin(p), end(p));\n       \
    \ }\n\n        std::vector<T> eval(std::vector<poly> &tree, int v, auto l, auto\
    \ r) { // auxiliary evaluation function\n            if(r - l == 1) {\n      \
    \          return {eval(*l)};\n            } else {\n                auto m =\
    \ l + (r - l) / 2;\n                auto A = (*this % tree[2 * v]).eval(tree,\
    \ 2 * v, l, m);\n                auto B = (*this % tree[2 * v + 1]).eval(tree,\
    \ 2 * v + 1, m, r);\n                A.insert(end(A), begin(B), end(B));\n   \
    \             return A;\n            }\n        }\n        \n        std::vector<T>\
    \ eval(std::vector<T> x) { // evaluate polynomial in (x1, ..., xn)\n         \
    \   int n = x.size();\n            if(is_zero()) {\n                return std::vector<T>(n,\
    \ T(0));\n            }\n            std::vector<poly> tree(4 * n);\n        \
    \    build(tree, 1, begin(x), end(x));\n            return eval(tree, 1, begin(x),\
    \ end(x));\n        }\n        \n        poly inter(std::vector<poly> &tree, int\
    \ v, auto ly, auto ry) { // auxiliary interpolation function\n            if(ry\
    \ - ly == 1) {\n                return {*ly / a[0]};\n            } else {\n \
    \               auto my = ly + (ry - ly) / 2;\n                auto A = (*this\
    \ % tree[2 * v]).inter(tree, 2 * v, ly, my);\n                auto B = (*this\
    \ % tree[2 * v + 1]).inter(tree, 2 * v + 1, my, ry);\n                return A\
    \ * tree[2 * v + 1] + B * tree[2 * v];\n            }\n        }\n        \n \
    \       static auto inter(std::vector<T> x, std::vector<T> y) { // interpolates\
    \ minimum polynomial from (xi, yi) pairs\n            int n = x.size();\n    \
    \        std::vector<poly> tree(4 * n);\n            return build(tree, 1, begin(x),\
    \ end(x)).deriv().inter(tree, 1, begin(y), end(y));\n        }\n\n        static\
    \ auto resultant(poly a, poly b) { // computes resultant of a and b\n        \
    \    if(b.is_zero()) {\n                return 0;\n            } else if(b.deg()\
    \ == 0) {\n                return bpow(b.lead(), a.deg());\n            } else\
    \ {\n                int pw = a.deg();\n                a %= b;\n            \
    \    pw -= a.deg();\n                auto mul = bpow(b.lead(), pw) * T((b.deg()\
    \ & a.deg() & 1) ? -1 : 1);\n                auto ans = resultant(b, a);\n   \
    \             return ans * mul;\n            }\n        }\n                \n\
    \        static poly xk(size_t n) { // P(x) = x^n\n            return poly(T(1)).mul_xk(n);\n\
    \        }\n        \n        static poly ones(size_t n) { // P(x) = 1 + x + ...\
    \ + x^{n-1} \n            return std::vector<T>(n, 1);\n        }\n        \n\
    \        static poly expx(size_t n) { // P(x) = e^x (mod x^n)\n            return\
    \ ones(n).borel();\n        }\n\n        static poly log1px(size_t n) { // P(x)\
    \ = log(1+x) (mod x^n)\n            std::vector<T> coeffs(n, 0);\n           \
    \ for(size_t i = 1; i < n; i++) {\n                coeffs[i] = (i & 1 ? T(i).inv()\
    \ : -T(i).inv());\n            }\n            return coeffs;\n        }\n\n  \
    \      static poly log1mx(size_t n) { // P(x) = log(1-x) (mod x^n)\n         \
    \   return -ones(n).integr();\n        }\n        \n        // [x^k] (a corr b)\
    \ = sum_{i} a{(k-m)+i}*bi\n        static poly corr(poly a, poly b) { // cross-correlation\n\
    \            return a * b.reverse();\n        }\n\n        // [x^k] (a semicorr\
    \ b) = sum_i a{i+k} * b{i}\n        static poly semicorr(poly a, poly b) {\n \
    \           return corr(a, b).div_xk(b.deg());\n        }\n        \n        poly\
    \ invborel() const { // ak *= k!\n            auto res = *this;\n            for(int\
    \ i = 0; i <= deg(); i++) {\n                res.coef(i) *= fact<T>(i);\n    \
    \        }\n            return res;\n        }\n        \n        poly borel()\
    \ const { // ak /= k!\n            auto res = *this;\n            for(int i =\
    \ 0; i <= deg(); i++) {\n                res.coef(i) *= rfact<T>(i);\n       \
    \     }\n            return res;\n        }\n        \n        poly shift(T a)\
    \ const { // P(x + a)\n            return semicorr(invborel(), expx(deg() + 1).mulx(a)).borel();\n\
    \        }\n        \n        poly x2() { // P(x) -> P(x^2)\n            std::vector<T>\
    \ res(2 * a.size());\n            for(size_t i = 0; i < a.size(); i++) {\n   \
    \             res[2 * i] = a[i];\n            }\n            return res;\n   \
    \     }\n        \n        // Return {P0, P1}, where P(x) = P0(x) + xP1(x)\n \
    \       std::pair<poly, poly> bisect() const {\n            std::vector<T> res[2];\n\
    \            res[0].reserve(deg() / 2 + 1);\n            res[1].reserve(deg()\
    \ / 2 + 1);\n            for(int i = 0; i <= deg(); i++) {\n                res[i\
    \ % 2].push_back(a[i]);\n            }\n            return {res[0], res[1]};\n\
    \        }\n        \n        // Find [x^k] P / Q\n        static T kth_rec(poly\
    \ P, poly Q, int64_t k) {\n            while(k > Q.deg()) {\n                int\
    \ n = Q.a.size();\n                auto [Q0, Q1] = Q.mulx(-1).bisect();\n    \
    \            auto [P0, P1] = P.bisect();\n                \n                int\
    \ N = fft::com_size((n + 1) / 2, (n + 1) / 2);\n                \n           \
    \     auto Q0f = fft::dft(Q0.a, N);\n                auto Q1f = fft::dft(Q1.a,\
    \ N);\n                auto P0f = fft::dft(P0.a, N);\n                auto P1f\
    \ = fft::dft(P1.a, N);\n                \n                if(k % 2) {\n      \
    \              P = poly(Q0f * P1f) + poly(Q1f * P0f);\n                } else\
    \ {\n                    P = poly(Q0f * P0f) + poly(Q1f * P1f).mul_xk(1);\n  \
    \              }\n                Q = poly(Q0f * Q0f) - poly(Q1f * Q1f).mul_xk(1);\n\
    \                k /= 2;\n            }\n            return (P * Q.inv(Q.deg()\
    \ + 1))[k];\n        }\n        \n        poly inv(int n) const { // get inverse\
    \ series mod x^n\n            auto Q = mod_xk(n);\n            if(n == 1) {\n\
    \                return Q[0].inv();\n            }\n            // Q(-x) = P0(x^2)\
    \ + xP1(x^2)\n            auto [P0, P1] = Q.mulx(-1).bisect();\n            \n\
    \            int N = fft::com_size((n + 1) / 2, (n + 1) / 2);\n            \n\
    \            auto P0f = fft::dft(P0.a, N);\n            auto P1f = fft::dft(P1.a,\
    \ N);\n            \n            auto TTf = fft::dft(( // Q(x)*Q(-x) = Q0(x^2)^2\
    \ - x^2 Q1(x^2)^2\n                poly(P0f * P0f) - poly(P1f * P1f).mul_xk(1)\n\
    \            ).inv((n + 1) / 2).a, N);\n            \n            return (\n \
    \               poly(P0f * TTf).x2() + poly(P1f * TTf).x2().mul_xk(1)\n      \
    \      ).mod_xk(n);\n        }\n        \n        // compute A(B(x)) mod x^n in\
    \ O(n^2)\n        static poly compose(poly A, poly B, int n) {\n            int\
    \ q = std::sqrt(n);\n            std::vector<poly> Bk(q);\n            auto Bq\
    \ = B.pow(q, n);\n            Bk[0] = poly(T(1));\n            for(int i = 1;\
    \ i < q; i++) {\n                Bk[i] = (Bk[i - 1] * B).mod_xk(n);\n        \
    \    }\n            poly Bqk(1);\n            poly ans;\n            for(int i\
    \ = 0; i <= n / q; i++) {\n                poly cur;\n                for(int\
    \ j = 0; j < q; j++) {\n                    cur += Bk[j] * A[i * q + j];\n   \
    \             }\n                ans += (Bqk * cur).mod_xk(n);\n             \
    \   Bqk = (Bqk * Bq).mod_xk(n);\n            }\n            return ans;\n    \
    \    }\n        \n        // compute A(B(x)) mod x^n in O(sqrt(pqn log^3 n))\n\
    \        // preferrable when p = deg A and q = deg B\n        // are much less\
    \ than n\n        static poly compose_large(poly A, poly B, int n) {\n       \
    \     if(B[0] != T(0)) {\n                return compose_large(A.shift(B[0]),\
    \ B - B[0], n);\n            }\n            \n            int q = std::sqrt(n);\n\
    \            auto [B0, B1] = make_pair(B.mod_xk(q), B.div_xk(q));\n          \
    \  \n            B0 = B0.div_xk(1);\n            std::vector<poly> pw(A.deg()\
    \ + 1);\n            auto getpow = [&](int k) {\n                return pw[k].is_zero()\
    \ ? pw[k] = B0.pow(k, n - k) : pw[k];\n            };\n            \n        \
    \    std::function<poly(poly const&, int, int)> compose_dac = [&getpow, &compose_dac](poly\
    \ const& f, int m, int N) {\n                if(f.deg() <= 0) {\n            \
    \        return f;\n                }\n                int k = m / 2;\n      \
    \          auto [f0, f1] = make_pair(f.mod_xk(k), f.div_xk(k));\n            \
    \    auto [A, B] = make_pair(compose_dac(f0, k, N), compose_dac(f1, m - k, N -\
    \ k));\n                return (A + (B.mod_xk(N - k) * getpow(k).mod_xk(N - k)).mul_xk(k)).mod_xk(N);\n\
    \            };\n            \n            int r = n / q;\n            auto Ar\
    \ = A.deriv(r);\n            auto AB0 = compose_dac(Ar, Ar.deg() + 1, n);\n  \
    \          \n            auto Bd = B0.mul_xk(1).deriv();\n            \n     \
    \       poly ans = T(0);\n            \n            std::vector<poly> B1p(r +\
    \ 1);\n            B1p[0] = poly(T(1));\n            for(int i = 1; i <= r; i++)\
    \ {\n                B1p[i] = (B1p[i - 1] * B1.mod_xk(n - i * q)).mod_xk(n - i\
    \ * q);\n            }\n            while(r >= 0) {\n                ans += (AB0.mod_xk(n\
    \ - r * q) * rfact<T>(r) * B1p[r]).mul_xk(r * q).mod_xk(n);\n                r--;\n\
    \                if(r >= 0) {\n                    AB0 = ((AB0 * Bd).integr()\
    \ + A[r] * fact<T>(r)).mod_xk(n);\n                }\n            }\n        \
    \    \n            return ans;\n        }\n    };\n    \n    static auto operator\
    \ * (const auto& a, const poly<auto>& b) {\n        return b * a;\n    }\n};\n\
    \n#line 5 \"verify/algebra/polynomial/poly_sqrt.test.cpp\"\n#include <bits/stdc++.h>\n\
    \nusing namespace std;\nusing namespace cp_algo::algebra;\n\nconst int mod = 998244353;\n\
    typedef modular<mod> base;\ntypedef poly<base> polyn;\n\nvoid solve() {\n    int\
    \ n;\n    cin >> n;\n    vector<base> a(n);\n    copy_n(istream_iterator<base>(cin),\
    \ n, begin(a));\n    auto res = polyn(a).sqrt(n);\n    if(res) {\n        res->print(n);\n\
    \    } else {\n        cout << -1 << \"\\n\";\n    }\n}\n\nsigned main() {\n \
    \   //freopen(\"input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n  \
    \  cin.tie(0);\n    int t = 1;\n    while(t--) {\n        solve();\n    }\n}\n"
  code: "// @brief Sqrt of Formal Power Series\n#define PROBLEM \"https://judge.yosupo.jp/problem/sqrt_of_formal_power_series\"\
    \n#include \"cp-algo/algebra/polynomial.hpp\"\n#include \"cp-algo/algebra/modular.hpp\"\
    \n#include <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::algebra;\n\
    \nconst int mod = 998244353;\ntypedef modular<mod> base;\ntypedef poly<base> polyn;\n\
    \nvoid solve() {\n    int n;\n    cin >> n;\n    vector<base> a(n);\n    copy_n(istream_iterator<base>(cin),\
    \ n, begin(a));\n    auto res = polyn(a).sqrt(n);\n    if(res) {\n        res->print(n);\n\
    \    } else {\n        cout << -1 << \"\\n\";\n    }\n}\n\nsigned main() {\n \
    \   //freopen(\"input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n  \
    \  cin.tie(0);\n    int t = 1;\n    while(t--) {\n        solve();\n    }\n}\n"
  dependsOn:
  - cp-algo/algebra/polynomial.hpp
  - cp-algo/algebra/common.hpp
  - cp-algo/algebra/fft.hpp
  - cp-algo/algebra/modular.hpp
  - cp-algo/random/rng.hpp
  - cp-algo/algebra/affine.hpp
  - cp-algo/algebra/modular.hpp
  isVerificationFile: true
  path: verify/algebra/polynomial/poly_sqrt.test.cpp
  requiredBy: []
  timestamp: '2024-02-11 00:07:44+01:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/algebra/polynomial/poly_sqrt.test.cpp
layout: document
redirect_from:
- /verify/verify/algebra/polynomial/poly_sqrt.test.cpp
- /verify/verify/algebra/polynomial/poly_sqrt.test.cpp.html
title: Sqrt of Formal Power Series
---
