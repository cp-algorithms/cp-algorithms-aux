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
  - icon: ':heavy_check_mark:'
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
    path: verify/algebra/polynomial/find_linrec.test.cpp
    title: Find Linear Recurrence
  - icon: ':heavy_check_mark:'
    path: verify/algebra/polynomial/poly_sqrt.test.cpp
    title: Sqrt of Formal Power Series
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  - icon: ':heavy_check_mark:'
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
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n\
    #include <cassert>\nnamespace cp_algo::algebra {\n    // a * x + b\n    template<typename\
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
    \ * A + b * B, c * A + d * B};\n        }\n    };\n}\n\n"
  code: "#ifndef CP_ALGO_ALGEBRA_AFFINE_HPP\n#define CP_ALGO_ALGEBRA_AFFINE_HPP\n\
    #include <optional>\n#include <cassert>\nnamespace cp_algo::algebra {\n    //\
    \ a * x + b\n    template<typename base>\n    struct lin {\n        base a = 1,\
    \ b = 0;\n        std::optional<base> c;\n        lin() {}\n        lin(base b):\
    \ a(0), b(b) {}\n        lin(base a, base b): a(a), b(b) {}\n        lin(base\
    \ a, base b, base _c): a(a), b(b), c(_c) {}\n\n        // polynomial product modulo\
    \ x^2 - c\n        lin operator * (const lin& t) {\n            assert(c && t.c\
    \ && *c == *t.c);\n            return {a * t.b + b * t.a, b * t.b + a * t.a *\
    \ (*c), *c};\n        }\n\n        // a * (t.a * x + t.b) + b\n        lin apply(lin\
    \ const& t) const {\n            return {a * t.a, a * t.b + b};\n        }\n\n\
    \        void prepend(lin const& t) {\n            *this = t.apply(*this);\n \
    \       }\n\n        base eval(base x) const {\n            return a * x + b;\n\
    \        }\n    };\n\n    // (ax+b) / (cx+d)\n    template<typename base>\n  \
    \  struct linfrac {\n        // default constructor for a continued fraction block\n\
    \        base a, b = base(1), c = base(1), d = base(0);\n        linfrac(base\
    \ a): a(a) {}\n        linfrac(base a, base b, base c, base d): a(a), b(b), c(c),\
    \ d(d) {}\n        \n        // composition of two linfracs\n        linfrac operator\
    \ *(linfrac const& t) {\n            auto [A, C] = apply(t.a, t.c);\n        \
    \    auto [B, D] = apply(t.b, t.d);\n            return {A, B, C, D};\n      \
    \  }\n        \n        linfrac adj() {\n            return {d, -b, -c, a};\n\
    \        }\n        \n        // apply linfrac to A/B\n        auto apply(base\
    \ A, base B) {\n            return std::pair{a * A + b * B, c * A + d * B};\n\
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
  timestamp: '2024-02-11 12:35:24+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
  - verify/data_structures/treap/range_reverse_range_sum.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  - verify/algebra/matrix/matrix_pow.test.cpp
  - verify/algebra/polynomial/convolution107.test.cpp
  - verify/algebra/polynomial/find_linrec.test.cpp
  - verify/algebra/polynomial/poly_sqrt.test.cpp
documentation_of: cp-algo/algebra/affine.hpp
layout: document
redirect_from:
- /library/cp-algo/algebra/affine.hpp
- /library/cp-algo/algebra/affine.hpp.html
title: cp-algo/algebra/affine.hpp
---
