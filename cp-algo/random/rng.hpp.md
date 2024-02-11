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
    path: cp-algo/data_structures/treap.hpp
    title: cp-algo/data_structures/treap.hpp
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
    path: verify/data_structures/treap/cartesian_tree.test.cpp
    title: Build Cartesian Tree
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
  bundledCode: "#line 1 \"cp-algo/random/rng.hpp\"\n\n\n#include <chrono>\n#include\
    \ <random>\nnamespace cp_algo::random {\n    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());\
    \ \n}\n\n"
  code: "#ifndef CP_ALGO_RANDOM_RNG_HPP\n#define CP_ALGO_RANDOM_RNG_HPP\n#include\
    \ <chrono>\n#include <random>\nnamespace cp_algo::random {\n    std::mt19937_64\
    \ rng(std::chrono::steady_clock::now().time_since_epoch().count()); \n}\n#endif\
    \ // CP_ALGO_RANDOM_RNG_HPP"
  dependsOn: []
  isVerificationFile: false
  path: cp-algo/random/rng.hpp
  requiredBy:
  - cp-algo/data_structures/treap.hpp
  - cp-algo/algebra/matrix.hpp
  - cp-algo/algebra/polynomial.hpp
  - cp-algo/algebra/modular.hpp
  - cp-algo/algebra/fft.hpp
  timestamp: '2024-02-11 12:35:24+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
  - verify/data_structures/treap/range_reverse_range_sum.test.cpp
  - verify/data_structures/treap/cartesian_tree.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  - verify/algebra/matrix/matrix_pow.test.cpp
  - verify/algebra/polynomial/convolution107.test.cpp
  - verify/algebra/polynomial/find_linrec.test.cpp
  - verify/algebra/polynomial/poly_sqrt.test.cpp
documentation_of: cp-algo/random/rng.hpp
layout: document
redirect_from:
- /library/cp-algo/random/rng.hpp
- /library/cp-algo/random/rng.hpp.html
title: cp-algo/random/rng.hpp
---
