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
  bundledCode: "#line 1 \"src/algebra/common.cpp\"\nnamespace algebra { // common\n\
    \    const int maxn = 1 << 20;\n    const int magic = 250; // threshold for sizes\
    \ to run the naive algo\n    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());\
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
    \  }\n}\n"
  code: "namespace algebra { // common\n    const int maxn = 1 << 20;\n    const int\
    \ magic = 250; // threshold for sizes to run the naive algo\n    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());\
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
    \  }\n}\n"
  dependsOn: []
  isVerificationFile: false
  path: src/algebra/common.cpp
  requiredBy: []
  timestamp: '2024-02-10 16:40:11+01:00'
  verificationStatus: LIBRARY_NO_TESTS
  verifiedWith: []
documentation_of: src/algebra/common.cpp
layout: document
redirect_from:
- /library/src/algebra/common.cpp
- /library/src/algebra/common.cpp.html
title: src/algebra/common.cpp
---
