---
data:
  _extendedDependsOn: []
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/fft.hpp
    title: cp-algo/algebra/fft.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/matrix.hpp
    title: cp-algo/algebra/matrix.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/polynomial.hpp
    title: cp-algo/algebra/polynomial.hpp
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/algebra/convolution107.test.cpp
    title: Convolution mod $10^9+7$
  - icon: ':heavy_check_mark:'
    path: verify/algebra/matrix_pow.test.cpp
    title: Pow of Matrix
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/algebra/common.hpp\"\n\n\n#include <cstdint>\nnamespace\
    \ cp_algo::algebra {\n    const int maxn = 1 << 20;\n    const int magic = 250;\
    \ // threshold for sizes to run the naive algo\n\n    auto bpow(auto x, int64_t\
    \ n, auto ans) {\n        for(; n; n /= 2, x = x * x) {\n            if(n % 2)\
    \ {\n                ans = ans * x;\n            }\n        }\n        return\
    \ ans;\n    }\n    template<typename T>\n    T bpow(T const& x, int64_t n) {\n\
    \        return bpow(x, n, T(1));\n    }\n\n    template<typename T>\n    T fact(int\
    \ n) {\n        static T F[maxn];\n        static bool init = false;\n       \
    \ if(!init) {\n            F[0] = T(1);\n            for(int i = 1; i < maxn;\
    \ i++) {\n                F[i] = F[i - 1] * T(i);\n            }\n           \
    \ init = true;\n        }\n        return F[n];\n    }\n    \n    template<typename\
    \ T>\n    T rfact(int n) {\n        static T F[maxn];\n        static bool init\
    \ = false;\n        if(!init) {\n            F[maxn - 1] = T(1) / fact<T>(maxn\
    \ - 1);\n            for(int i = maxn - 2; i >= 0; i--) {\n                F[i]\
    \ = F[i + 1] * T(i + 1);\n            }\n            init = true;\n        }\n\
    \        return F[n];\n    }\n\n    template<typename T>\n    T small_inv(int\
    \ n) {\n        static T F[maxn];\n        static bool init = false;\n       \
    \ if(!init) {\n            for(int i = 1; i < maxn; i++) {\n                F[i]\
    \ = rfact<T>(i) * fact<T>(i - 1);\n            }\n            init = true;\n \
    \       }\n        return F[n];\n    }\n}\n\n"
  code: "#ifndef CP_ALGO_ALGEBRA_COMMON_HPP\n#define CP_ALGO_ALGEBRA_COMMON_HPP\n\
    #include <cstdint>\nnamespace cp_algo::algebra {\n    const int maxn = 1 << 20;\n\
    \    const int magic = 250; // threshold for sizes to run the naive algo\n\n \
    \   auto bpow(auto x, int64_t n, auto ans) {\n        for(; n; n /= 2, x = x *\
    \ x) {\n            if(n % 2) {\n                ans = ans * x;\n            }\n\
    \        }\n        return ans;\n    }\n    template<typename T>\n    T bpow(T\
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
    \     init = true;\n        }\n        return F[n];\n    }\n}\n#endif // CP_ALGO_ALGEBRA_COMMON_HPP\n"
  dependsOn: []
  isVerificationFile: false
  path: cp-algo/algebra/common.hpp
  requiredBy:
  - cp-algo/algebra/matrix.hpp
  - cp-algo/algebra/polynomial.hpp
  - cp-algo/algebra/modular.hpp
  - cp-algo/algebra/fft.hpp
  timestamp: '2024-02-10 22:44:24+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
  - verify/algebra/matrix_pow.test.cpp
  - verify/algebra/convolution107.test.cpp
documentation_of: cp-algo/algebra/common.hpp
layout: document
redirect_from:
- /library/cp-algo/algebra/common.hpp
- /library/cp-algo/algebra/common.hpp.html
title: cp-algo/algebra/common.hpp
---
