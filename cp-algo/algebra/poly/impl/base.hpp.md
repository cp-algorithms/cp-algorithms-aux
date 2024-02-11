---
data:
  _extendedDependsOn: []
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
  bundledCode: "#line 1 \"cp-algo/algebra/poly/impl/base.hpp\"\n\n\n#include <functional>\n\
    #include <algorithm>\n// really basic operations, typically taking O(n)\nnamespace\
    \ cp_algo::algebra::poly::impl {\n    void normalize(auto& p) {\n        while(p.deg()\
    \ >= 0 && p.lead() == 0) {\n            p.a.pop_back();\n        }\n    }\n  \
    \  auto neg(auto const& p) {\n        auto res = p.a;\n        std::ranges::transform(res,\
    \ begin(res), std::negate<>{});\n        return res;\n    }\n}\n\n"
  code: "#ifndef CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP\n#define CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP\n\
    #include <functional>\n#include <algorithm>\n// really basic operations, typically\
    \ taking O(n)\nnamespace cp_algo::algebra::poly::impl {\n    void normalize(auto&\
    \ p) {\n        while(p.deg() >= 0 && p.lead() == 0) {\n            p.a.pop_back();\n\
    \        }\n    }\n    auto neg(auto const& p) {\n        auto res = p.a;\n  \
    \      std::ranges::transform(res, begin(res), std::negate<>{});\n        return\
    \ res;\n    }\n}\n#endif // CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP\n"
  dependsOn: []
  isVerificationFile: false
  path: cp-algo/algebra/poly/impl/base.hpp
  requiredBy:
  - cp-algo/algebra/poly.hpp
  timestamp: '2024-02-11 14:42:51+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/algebra/poly/convolution107.test.cpp
  - verify/algebra/poly/find_linrec.test.cpp
  - verify/algebra/poly/poly_sqrt.test.cpp
documentation_of: cp-algo/algebra/poly/impl/base.hpp
layout: document
redirect_from:
- /library/cp-algo/algebra/poly/impl/base.hpp
- /library/cp-algo/algebra/poly/impl/base.hpp.html
title: cp-algo/algebra/poly/impl/base.hpp
---
