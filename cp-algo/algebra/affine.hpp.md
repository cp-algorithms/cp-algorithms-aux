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
  - icon: ':question:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/polynomial.hpp
    title: cp-algo/algebra/polynomial.hpp
  - icon: ':x:'
    path: cp-algo/data_structures/segment_tree/metas/affine.hpp
    title: cp-algo/data_structures/segment_tree/metas/affine.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/metas/reverse.hpp
    title: cp-algo/data_structures/treap/metas/reverse.hpp
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/algebra/matrix/matrix_pow.test.cpp
    title: Pow of Matrix
  - icon: ':heavy_check_mark:'
    path: verify/algebra/polynomial/convolution107.test.cpp
    title: Convolution mod $10^9+7$
  - icon: ':heavy_check_mark:'
    path: verify/algebra/polynomial/poly_sqrt.test.cpp
    title: Sqrt of Formal Power Series
  - icon: ':x:'
    path: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  - icon: ':x:'
    path: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
    title: Dynamic Range Affine Range Sum
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
    title: Dynamic Range Affine Range Sum
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/treap/range_reverse_range_sum.test.cpp
    title: Range Reverse Range Sum
  _isVerificationFailed: true
  _pathExtension: hpp
  _verificationStatusIcon: ':question:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n\
    #include <cassert>\nnamespace cp_algo::algebra {\n    template<typename base>\n\
    \    // a * x + b\n    struct lin {\n        base a = 1, b = 0;\n        std::optional<base>\
    \ c;\n        lin() {}\n        lin(base b): a(0), b(b) {}\n        lin(base a,\
    \ base b): a(a), b(b) {}\n        lin(base a, base b, base _c): a(a), b(b), c(_c)\
    \ {}\n\n        // polynomial product modulo x^2 - c\n        lin operator * (const\
    \ lin& t) {\n            assert(c && t.c && *c == *t.c);\n            return lin(a\
    \ * t.b + b * t.a, b * t.b + a * t.a * (*c), *c);\n        }\n\n        // a *\
    \ (t.a * x + t.b) + b\n        lin compose(lin const& t) const {\n           \
    \ return lin{a * t.a, a * t.b + b};\n        }\n\n        void prepend(lin const&\
    \ t) {\n            *this = t.compose(*this);\n        }\n\n        base eval(base\
    \ x) const {\n            return a * x + b;\n        }\n    };\n}\n\n"
  code: "#ifndef CP_ALGO_ALGEBRA_AFFINE_HPP\n#define CP_ALGO_ALGEBRA_AFFINE_HPP\n\
    #include <optional>\n#include <cassert>\nnamespace cp_algo::algebra {\n    template<typename\
    \ base>\n    // a * x + b\n    struct lin {\n        base a = 1, b = 0;\n    \
    \    std::optional<base> c;\n        lin() {}\n        lin(base b): a(0), b(b)\
    \ {}\n        lin(base a, base b): a(a), b(b) {}\n        lin(base a, base b,\
    \ base _c): a(a), b(b), c(_c) {}\n\n        // polynomial product modulo x^2 -\
    \ c\n        lin operator * (const lin& t) {\n            assert(c && t.c && *c\
    \ == *t.c);\n            return lin(a * t.b + b * t.a, b * t.b + a * t.a * (*c),\
    \ *c);\n        }\n\n        // a * (t.a * x + t.b) + b\n        lin compose(lin\
    \ const& t) const {\n            return lin{a * t.a, a * t.b + b};\n        }\n\
    \n        void prepend(lin const& t) {\n            *this = t.compose(*this);\n\
    \        }\n\n        base eval(base x) const {\n            return a * x + b;\n\
    \        }\n    };\n}\n#endif // CP_ALGO_ALGEBRA_AFFINE_HPP"
  dependsOn: []
  isVerificationFile: false
  path: cp-algo/algebra/affine.hpp
  requiredBy:
  - cp-algo/data_structures/segment_tree/metas/affine.hpp
  - cp-algo/data_structures/treap/metas/reverse.hpp
  - cp-algo/algebra/matrix.hpp
  - cp-algo/algebra/polynomial.hpp
  - cp-algo/algebra/modular.hpp
  - cp-algo/algebra/fft.hpp
  timestamp: '2024-02-10 23:55:00+01:00'
  verificationStatus: LIBRARY_SOME_WA
  verifiedWith:
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
  - verify/data_structures/treap/range_reverse_range_sum.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  - verify/algebra/matrix/matrix_pow.test.cpp
  - verify/algebra/polynomial/convolution107.test.cpp
  - verify/algebra/polynomial/poly_sqrt.test.cpp
documentation_of: cp-algo/algebra/affine.hpp
layout: document
redirect_from:
- /library/cp-algo/algebra/affine.hpp
- /library/cp-algo/algebra/affine.hpp.html
title: cp-algo/algebra/affine.hpp
---
