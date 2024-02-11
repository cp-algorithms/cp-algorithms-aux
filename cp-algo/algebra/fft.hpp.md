---
data:
  _extendedDependsOn:
  - icon: ':question:'
    path: cp-algo/algebra/affine.hpp
    title: cp-algo/algebra/affine.hpp
  - icon: ':question:'
    path: cp-algo/algebra/common.hpp
    title: cp-algo/algebra/common.hpp
  - icon: ':question:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
  - icon: ':question:'
    path: cp-algo/random/rng.hpp
    title: cp-algo/random/rng.hpp
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/poly.hpp
    title: cp-algo/algebra/poly.hpp
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/algebra/poly/convolution107.test.cpp
    title: Convolution mod $10^9+7$
  - icon: ':heavy_check_mark:'
    path: verify/algebra/poly/find_linrec.test.cpp
    title: Find Linear Recurrence
  - icon: ':heavy_check_mark:'
    path: verify/algebra/poly/poly_sqrt.test.cpp
    title: Sqrt of Formal Power Series
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/algebra/fft.hpp\"\n\n\n#line 1 \"cp-algo/algebra/common.hpp\"\
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
    \ = A * dft<m>(b, n);\n        }\n    }\n}\n\n"
  code: "#ifndef CP_ALGO_ALGEBRA_FFT_HPP\n#define CP_ALGO_ALGEBRA_FFT_HPP\n#include\
    \ \"common.hpp\"\n#include \"modular.hpp\"\n#include <algorithm>\n#include <cassert>\n\
    #include <vector>\nnamespace cp_algo::algebra::fft {\n    using ftype = double;\n\
    \    struct point {\n        ftype x, y;\n        \n        ftype real() {return\
    \ x;}\n        ftype imag() {return y;}\n        \n        point(): x(0), y(0){}\n\
    \        point(ftype x, ftype y = 0): x(x), y(y){}\n        \n        static point\
    \ polar(ftype rho, ftype ang) {\n            return point{rho * cos(ang), rho\
    \ * sin(ang)};\n        }\n        \n        point conj() const {\n          \
    \  return {x, -y};\n        }\n        \n        point operator +=(const point\
    \ &t) {x += t.x, y += t.y; return *this;}\n        point operator +(const point\
    \ &t) const {return point(*this) += t;}\n        point operator -(const point\
    \ &t) const {return {x - t.x, y - t.y};}\n        point operator *(const point\
    \ &t) const {return {x * t.x - y * t.y, x * t.y + y * t.x};}\n    };\n\n    point\
    \ w[maxn]; // w[2^n + k] = exp(pi * k / (2^n))\n    int bitr[maxn];// b[2^n +\
    \ k] = bitreverse(k)\n    const ftype pi = acos(-1);\n    bool initiated = 0;\n\
    \    void init() {\n        if(!initiated) {\n            for(int i = 1; i < maxn;\
    \ i *= 2) {\n                int ti = i / 2;\n                for(int j = 0; j\
    \ < i; j++) {\n                    w[i + j] = point::polar(ftype(1), pi * j /\
    \ i);\n                    if(ti) {\n                        bitr[i + j] = 2 *\
    \ bitr[ti + j % ti] + (j >= ti);\n                    }\n                }\n \
    \           }\n            initiated = 1;\n        }\n    }\n    \n    void fft(auto\
    \ &a, int n) {\n        init();\n        if(n == 1) {\n            return;\n \
    \       }\n        int hn = n / 2;\n        for(int i = 0; i < n; i++) {\n   \
    \         int ti = 2 * bitr[hn + i % hn] + (i > hn);\n            if(i < ti) {\n\
    \                std::swap(a[i], a[ti]);\n            }\n        }\n        for(int\
    \ i = 1; i < n; i *= 2) {\n            for(int j = 0; j < n; j += 2 * i) {\n \
    \               for(int k = j; k < j + i; k++) {\n                    point t\
    \ = a[k + i] * w[i + k - j];\n                    a[k + i] = a[k] - t;\n     \
    \               a[k] += t;\n                }\n            }\n        }\n    }\n\
    \    \n    void mul_slow(std::vector<auto> &a, const std::vector<auto> &b) {\n\
    \        if(a.empty() || b.empty()) {\n            a.clear();\n        } else\
    \ {\n            int n = a.size();\n            int m = b.size();\n          \
    \  a.resize(n + m - 1);\n            for(int k = n + m - 2; k >= 0; k--) {\n \
    \               a[k] *= b[0];\n                for(int j = std::max(k - n + 1,\
    \ 1); j < std::min(k + 1, m); j++) {\n                    a[k] += a[k - j] * b[j];\n\
    \                }\n            }\n        }\n    }\n    \n    template<int m>\n\
    \    struct dft {\n        static constexpr int split = 1 << 15;\n        std::vector<point>\
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
    \ = A * dft<m>(b, n);\n        }\n    }\n}\n#endif // CP_ALGO_ALGEBRA_FFT_HPP\n"
  dependsOn:
  - cp-algo/algebra/common.hpp
  - cp-algo/algebra/modular.hpp
  - cp-algo/random/rng.hpp
  - cp-algo/algebra/affine.hpp
  isVerificationFile: false
  path: cp-algo/algebra/fft.hpp
  requiredBy:
  - cp-algo/algebra/poly.hpp
  timestamp: '2024-02-11 14:42:51+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/algebra/poly/convolution107.test.cpp
  - verify/algebra/poly/find_linrec.test.cpp
  - verify/algebra/poly/poly_sqrt.test.cpp
documentation_of: cp-algo/algebra/fft.hpp
layout: document
redirect_from:
- /library/cp-algo/algebra/fft.hpp
- /library/cp-algo/algebra/fft.hpp.html
title: cp-algo/algebra/fft.hpp
---
