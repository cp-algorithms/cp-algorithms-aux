---
data:
  _extendedDependsOn: []
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segtree/metas/affine.hpp
    title: cp-algo/data_structures/segtree/metas/affine.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp
    title: cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/segtree/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/segtree/range_chmin_chmax_add_range_sum.test.cpp
    title: Range Chmin Chmax Add Range Sum
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/data_structures/segtree/metas/base.hpp\"\n\n\nnamespace\
    \ cp_algo::data_structures::segtree::metas {\n    template<typename derived_meta>\n\
    \    struct base_meta {\n        using meta = derived_meta;\n        virtual void\
    \ pull(meta const&, meta const&, int, int) {};\n        virtual void push(meta*,\
    \ meta*, int, int) {};\n    };\n}\n\n"
  code: "#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP\n#define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP\n\
    namespace cp_algo::data_structures::segtree::metas {\n    template<typename derived_meta>\n\
    \    struct base_meta {\n        using meta = derived_meta;\n        virtual void\
    \ pull(meta const&, meta const&, int, int) {};\n        virtual void push(meta*,\
    \ meta*, int, int) {};\n    };\n}\n#endif // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP"
  dependsOn: []
  isVerificationFile: false
  path: cp-algo/data_structures/segtree/metas/base.hpp
  requiredBy:
  - cp-algo/data_structures/segtree/metas/affine.hpp
  - cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp
  timestamp: '2024-02-11 15:34:32+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/segtree/range_affine_range_sum.test.cpp
  - verify/data_structures/segtree/range_chmin_chmax_add_range_sum.test.cpp
documentation_of: cp-algo/data_structures/segtree/metas/base.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/segtree/metas/base.hpp
- /library/cp-algo/data_structures/segtree/metas/base.hpp.html
title: cp-algo/data_structures/segtree/metas/base.hpp
---
