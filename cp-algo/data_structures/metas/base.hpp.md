---
data:
  _extendedDependsOn: []
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/metas/affine.hpp
    title: cp-algo/data_structures/metas/affine.hpp
  - icon: ':warning:'
    path: cp-algo/data_structures/metas/all.hpp
    title: cp-algo/data_structures/metas/all.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/metas/chmin_chmax_add.hpp
    title: cp-algo/data_structures/metas/chmin_chmax_add.hpp
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/range_chmin_chmax_add_range_sum.test.cpp
    title: Range Chmin Chmax Add Range Sum
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/data_structures/metas/base.hpp\"\n\n\n#include <functional>\n\
    #include <algorithm>\n#include <cstdint>\nnamespace data_structures {\n    namespace\
    \ segment_tree {\n        namespace metas {\n            template<typename derived_meta>\n\
    \            struct base_meta {\n                using meta = derived_meta;\n\
    \                virtual void pull(meta const&, meta const&, int, int) = 0;\n\
    \                virtual void push(meta*, meta*, int, int) = 0;\n            };\n\
    \        }\n    }\n}\n\n"
  code: "#ifndef DATA_STRUCTURES_METAS_BASE_HPP\n#define DATA_STRUCTURES_METAS_BASE_HPP\n\
    #include <functional>\n#include <algorithm>\n#include <cstdint>\nnamespace data_structures\
    \ {\n    namespace segment_tree {\n        namespace metas {\n            template<typename\
    \ derived_meta>\n            struct base_meta {\n                using meta =\
    \ derived_meta;\n                virtual void pull(meta const&, meta const&, int,\
    \ int) = 0;\n                virtual void push(meta*, meta*, int, int) = 0;\n\
    \            };\n        }\n    }\n}\n#endif // DATA_STRUCTURES_METAS_BASE_HPP"
  dependsOn: []
  isVerificationFile: false
  path: cp-algo/data_structures/metas/base.hpp
  requiredBy:
  - cp-algo/data_structures/metas/affine.hpp
  - cp-algo/data_structures/metas/all.hpp
  - cp-algo/data_structures/metas/chmin_chmax_add.hpp
  timestamp: '2024-02-10 20:45:15+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/range_affine_range_sum.test.cpp
  - verify/data_structures/range_chmin_chmax_add_range_sum.test.cpp
documentation_of: cp-algo/data_structures/metas/base.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/metas/base.hpp
- /library/cp-algo/data_structures/metas/base.hpp.html
title: cp-algo/data_structures/metas/base.hpp
---
