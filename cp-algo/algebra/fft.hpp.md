---
data:
  _extendedDependsOn:
  - icon: ':x:'
    path: cp-algo/algebra/common.hpp
    title: cp-algo/algebra/common.hpp
  - icon: ':x:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
  _extendedRequiredBy:
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
    links: []
  bundledCode: "Traceback (most recent call last):\n  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/documentation/build.py\"\
    , line 71, in _render_source_code_stat\n    bundled_code = language.bundle(stat.path,\
    \ basedir=basedir, options={'include_paths': [basedir]}).decode()\n  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus.py\"\
    , line 187, in bundle\n    bundler.update(path)\n  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 401, in update\n    self.update(self._resolve(pathlib.Path(included), included_from=path))\n\
    \  File \"/home/runner/.local/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 260, in _resolve\n    raise BundleErrorAt(path, -1, \"no such header\"\
    )\nonlinejudge_verify.languages.cplusplus_bundle.BundleErrorAt: cassert: line\
    \ -1: no such header\n"
  code: "#ifndef ALGEBRA_FFT_HPP\n#define ALGEBRA_FFT_HPP\n#include \"common.hpp\"\
    \n#include \"modular.hpp\"\n#include <algorithm>\n#include \"cassert\"\n#include\
    \ <vector>\nnamespace algebra {\n    namespace fft {\n        using ftype = double;\n\
    \        struct point {\n            ftype x, y;\n            \n            ftype\
    \ real() {return x;}\n            ftype imag() {return y;}\n            \n   \
    \         point(): x(0), y(0){}\n            point(ftype x, ftype y = 0): x(x),\
    \ y(y){}\n            \n            static point polar(ftype rho, ftype ang) {\n\
    \                return point{rho * cos(ang), rho * sin(ang)};\n            }\n\
    \            \n            point conj() const {\n                return {x, -y};\n\
    \            }\n            \n            point operator +=(const point &t) {x\
    \ += t.x, y += t.y; return *this;}\n            point operator +(const point &t)\
    \ const {return point(*this) += t;}\n            point operator -(const point\
    \ &t) const {return {x - t.x, y - t.y};}\n            point operator *(const point\
    \ &t) const {return {x * t.x - y * t.y, x * t.y + y * t.x};}\n        };\n\n \
    \       point w[maxn]; // w[2^n + k] = exp(pi * k / (2^n))\n        int bitr[maxn];//\
    \ b[2^n + k] = bitreverse(k)\n        const ftype pi = acos(-1);\n        bool\
    \ initiated = 0;\n        void init() {\n            if(!initiated) {\n      \
    \          for(int i = 1; i < maxn; i *= 2) {\n                    int ti = i\
    \ / 2;\n                    for(int j = 0; j < i; j++) {\n                   \
    \     w[i + j] = point::polar(ftype(1), pi * j / i);\n                       \
    \ if(ti) {\n                            bitr[i + j] = 2 * bitr[ti + j % ti] +\
    \ (j >= ti);\n                        }\n                    }\n             \
    \   }\n                initiated = 1;\n            }\n        }\n        \n  \
    \      void fft(auto &a, int n) {\n            init();\n            if(n == 1)\
    \ {\n                return;\n            }\n            int hn = n / 2;\n   \
    \         for(int i = 0; i < n; i++) {\n                int ti = 2 * bitr[hn +\
    \ i % hn] + (i > hn);\n                if(i < ti) {\n                    std::swap(a[i],\
    \ a[ti]);\n                }\n            }\n            for(int i = 1; i < n;\
    \ i *= 2) {\n                for(int j = 0; j < n; j += 2 * i) {\n           \
    \         for(int k = j; k < j + i; k++) {\n                        point t =\
    \ a[k + i] * w[i + k - j];\n                        a[k + i] = a[k] - t;\n   \
    \                     a[k] += t;\n                    }\n                }\n \
    \           }\n        }\n        \n        void mul_slow(std::vector<auto> &a,\
    \ const std::vector<auto> &b) {\n            if(a.empty() || b.empty()) {\n  \
    \              a.clear();\n            } else {\n                int n = a.size();\n\
    \                int m = b.size();\n                a.resize(n + m - 1);\n   \
    \             for(int k = n + m - 2; k >= 0; k--) {\n                    a[k]\
    \ *= b[0];\n                    for(int j = std::max(k - n + 1, 1); j < std::min(k\
    \ + 1, m); j++) {\n                        a[k] += a[k - j] * b[j];\n        \
    \            }\n                }\n            }\n        }\n        \n      \
    \  template<int m>\n        struct dft {\n            static constexpr int split\
    \ = 1 << 15;\n            std::vector<point> A;\n            \n            dft(std::vector<modular<m>>\
    \ const& a, size_t n): A(n) {\n                for(size_t i = 0; i < std::min(n,\
    \ a.size()); i++) {\n                    A[i] = point(\n                     \
    \   a[i].rem() % split,\n                        a[i].rem() / split\n        \
    \            );\n                }\n                if(n) {\n                \
    \    fft(A, n);\n                }\n            }\n        \n            auto\
    \ operator * (dft const& B) {\n                assert(A.size() == B.A.size());\n\
    \                size_t n = A.size();\n                if(!n) {\n            \
    \        return std::vector<modular<m>>();\n                }\n              \
    \  std::vector<point> C(n), D(n);\n                for(size_t i = 0; i < n; i++)\
    \ {\n                    C[i] = A[i] * (B[i] + B[(n - i) % n].conj());\n     \
    \               D[i] = A[i] * (B[i] - B[(n - i) % n].conj());\n              \
    \  }\n                fft(C, n);\n                fft(D, n);\n               \
    \ reverse(begin(C) + 1, end(C));\n                reverse(begin(D) + 1, end(D));\n\
    \                int t = 2 * n;\n                std::vector<modular<m>> res(n);\n\
    \                for(size_t i = 0; i < n; i++) {\n                    modular<m>\
    \ A0 = llround(C[i].real() / t);\n                    modular<m> A1 = llround(C[i].imag()\
    \ / t + D[i].imag() / t);\n                    modular<m> A2 = llround(D[i].real()\
    \ / t);\n                    res[i] = A0 + A1 * split - A2 * split * split;\n\
    \                }\n                return res;\n            }\n            \n\
    \            point& operator [](int i) {return A[i];}\n            point operator\
    \ [](int i) const {return A[i];}\n        };\n        \n        size_t com_size(size_t\
    \ as, size_t bs) {\n            if(!as || !bs) {\n                return 0;\n\
    \            }\n            size_t n = as + bs - 1;\n            while(__builtin_popcount(n)\
    \ != 1) {\n                n++;\n            }\n            return n;\n      \
    \  }\n        \n        template<int m>\n        void mul(std::vector<modular<m>>\
    \ &a, std::vector<modular<m>> b) {\n            if(std::min(a.size(), b.size())\
    \ < magic) {\n                mul_slow(a, b);\n                return;\n     \
    \       }\n            auto n = com_size(a.size(), b.size());\n            auto\
    \ A = dft<m>(a, n);\n            if(a == b) {\n                a = A * A;\n  \
    \          } else {\n                a = A * dft<m>(b, n);\n            }\n  \
    \      }\n    }\n}\n#endif // ALGEBRA_FFT_HPP\n"
  dependsOn:
  - cp-algo/algebra/common.hpp
  - cp-algo/algebra/modular.hpp
  isVerificationFile: false
  path: cp-algo/algebra/fft.hpp
  requiredBy:
  - cp-algo/algebra/polynomial.hpp
  timestamp: '2024-02-10 18:45:35+01:00'
  verificationStatus: LIBRARY_ALL_WA
  verifiedWith:
  - verify/algebra/convolution107.test.cpp
documentation_of: cp-algo/algebra/fft.hpp
layout: document
redirect_from:
- /library/cp-algo/algebra/fft.hpp
- /library/cp-algo/algebra/fft.hpp.html
title: cp-algo/algebra/fft.hpp
---
