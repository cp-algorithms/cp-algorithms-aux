---
data:
  _extendedDependsOn: []
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap.hpp
    title: cp-algo/data_structures/treap.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/metas/base.hpp
    title: cp-algo/data_structures/treap/metas/base.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/metas/reverse.hpp
    title: cp-algo/data_structures/treap/metas/reverse.hpp
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/treap/cartesian_tree.test.cpp
    title: Build Cartesian Tree
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
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/treap/range_reverse_range_sum.test.cpp
    title: Range Reverse Range Sum
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: '#line 1 "cp-algo/data_structures/treap/common.hpp"



    #define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())


    '
  code: '#ifndef CP_ALGO_DATA_STRUCTURES_TREAP_COMMON_HPP

    #define CP_ALGO_DATA_STRUCTURES_TREAP_COMMON_HPP

    #define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())

    #endif // CP_ALGO_DATA_STRUCTURES_TREAP_COMMON_HPP'
  dependsOn: []
  isVerificationFile: false
  path: cp-algo/data_structures/treap/common.hpp
  requiredBy:
  - cp-algo/data_structures/treap/metas/base.hpp
  - cp-algo/data_structures/treap/metas/reverse.hpp
  - cp-algo/data_structures/treap.hpp
  timestamp: '2024-02-11 00:07:44+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/treap/range_reverse_range_sum.test.cpp
  - verify/data_structures/treap/range_reverse_range_sum.test.cpp
  - verify/data_structures/treap/cartesian_tree.test.cpp
  - verify/data_structures/treap/cartesian_tree.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
documentation_of: cp-algo/data_structures/treap/common.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/treap/common.hpp
- /library/cp-algo/data_structures/treap/common.hpp.html
title: cp-algo/data_structures/treap/common.hpp
---
