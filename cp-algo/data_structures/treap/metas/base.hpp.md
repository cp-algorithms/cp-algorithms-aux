---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/common.hpp
    title: cp-algo/data_structures/treap/common.hpp
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/metas/reverse.hpp
    title: cp-algo/data_structures/treap/metas/reverse.hpp
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/treap/cartesian_tree.test.cpp
    title: Build Cartesian Tree
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
  bundledCode: "#line 1 \"cp-algo/data_structures/treap/metas/base.hpp\"\n\n\n#line\
    \ 1 \"cp-algo/data_structures/treap/common.hpp\"\n\n\n#define _safe(t, op) (t\
    \ ? t->op : typename std::remove_reference_t<decltype(t->op)>())\n\n#line 4 \"\
    cp-algo/data_structures/treap/metas/base.hpp\"\n#include <functional>\n#include\
    \ <algorithm>\n#include <cstdint>\n#define _safe_meta(i, op) _safe(i, _meta.op)\n\
    namespace cp_algo::data_structures::treap::metas {\n    struct base_meta {\n \
    \       void pull(auto const, auto const){}\n        void push(auto&, auto&){}\n\
    \    };\n}\n\n"
  code: "#ifndef CP_ALGO_DATA_STRUCTURES_TREAP_METAS_BASE_HPP\n#define CP_ALGO_DATA_STRUCTURES_TREAP_METAS_BASE_HPP\n\
    #include \"../common.hpp\"\n#include <functional>\n#include <algorithm>\n#include\
    \ <cstdint>\n#define _safe_meta(i, op) _safe(i, _meta.op)\nnamespace cp_algo::data_structures::treap::metas\
    \ {\n    struct base_meta {\n        void pull(auto const, auto const){}\n   \
    \     void push(auto&, auto&){}\n    };\n}\n#endif // CP_ALGO_DATA_STRUCTURES_TREAP_METAS_BASE_HPP"
  dependsOn:
  - cp-algo/data_structures/treap/common.hpp
  isVerificationFile: false
  path: cp-algo/data_structures/treap/metas/base.hpp
  requiredBy:
  - cp-algo/data_structures/treap/metas/reverse.hpp
  timestamp: '2024-02-11 12:35:24+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/treap/range_reverse_range_sum.test.cpp
  - verify/data_structures/treap/cartesian_tree.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
documentation_of: cp-algo/data_structures/treap/metas/base.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/treap/metas/base.hpp
- /library/cp-algo/data_structures/treap/metas/base.hpp.html
title: cp-algo/data_structures/treap/metas/base.hpp
---
