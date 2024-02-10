---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segment_tree/metas/base.hpp
    title: cp-algo/data_structures/segment_tree/metas/base.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/data_structures/segment_tree/metas/affine.hpp\"\n\
    \n\n#line 1 \"cp-algo/data_structures/segment_tree/metas/base.hpp\"\n\n\n#include\
    \ <functional>\n#include <algorithm>\n#include <cstdint>\nnamespace cp_algo::data_structures::segment_tree::metas\
    \ {\n    template<typename derived_meta>\n    struct base_meta {\n        using\
    \ meta = derived_meta;\n        virtual void pull(meta const&, meta const&, int,\
    \ int) = 0;\n        virtual void push(meta*, meta*, int, int) = 0;\n    };\n\
    }\n\n#line 4 \"cp-algo/data_structures/segment_tree/metas/affine.hpp\"\nnamespace\
    \ cp_algo::data_structures::segment_tree::metas {\n    template<typename base>\n\
    \    struct affine_meta: base_meta<affine_meta<base>> {\n        using meta =\
    \ affine_meta;\n        struct lin {\n            base a = 1, b = 0;\n       \
    \     lin() {}\n            lin(base a, base b): a(a), b(b){}\n\n            //\
    \ a * (t.a * x + t.b) + b\n            lin operator * (lin const& t) const {\n\
    \                return lin{a * t.a, a * t.b + b};\n            }\n\n        \
    \    base apply(base x) const {\n                return a * x + b;\n         \
    \   }\n        };\n\n        base sum = 0;\n        lin to_push = {};\n\n    \
    \    affine_meta() {}\n        affine_meta(base sum): sum(sum) {}\n\n        void\
    \ push(meta *L, meta *R, int l, int r) override {\n            if(to_push.a !=\
    \ 1 || to_push.b != 0) {\n                sum = to_push.a * sum + to_push.b *\
    \ (r - l);\n                if(r - l > 1) {\n                    L->to_push =\
    \ to_push * L->to_push;\n                    R->to_push = to_push * R->to_push;\n\
    \                }\n                to_push = {};\n            }\n        }\n\n\
    \        void pull(meta const& L, meta const& R, int, int) override {\n      \
    \      sum = L.sum + R.sum;\n        }\n    };\n}\n\n"
  code: "#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP\n#define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP\n\
    #include \"base.hpp\"\nnamespace cp_algo::data_structures::segment_tree::metas\
    \ {\n    template<typename base>\n    struct affine_meta: base_meta<affine_meta<base>>\
    \ {\n        using meta = affine_meta;\n        struct lin {\n            base\
    \ a = 1, b = 0;\n            lin() {}\n            lin(base a, base b): a(a),\
    \ b(b){}\n\n            // a * (t.a * x + t.b) + b\n            lin operator *\
    \ (lin const& t) const {\n                return lin{a * t.a, a * t.b + b};\n\
    \            }\n\n            base apply(base x) const {\n                return\
    \ a * x + b;\n            }\n        };\n\n        base sum = 0;\n        lin\
    \ to_push = {};\n\n        affine_meta() {}\n        affine_meta(base sum): sum(sum)\
    \ {}\n\n        void push(meta *L, meta *R, int l, int r) override {\n       \
    \     if(to_push.a != 1 || to_push.b != 0) {\n                sum = to_push.a\
    \ * sum + to_push.b * (r - l);\n                if(r - l > 1) {\n            \
    \        L->to_push = to_push * L->to_push;\n                    R->to_push =\
    \ to_push * R->to_push;\n                }\n                to_push = {};\n  \
    \          }\n        }\n\n        void pull(meta const& L, meta const& R, int,\
    \ int) override {\n            sum = L.sum + R.sum;\n        }\n    };\n}\n#endif\
    \ // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP"
  dependsOn:
  - cp-algo/data_structures/segment_tree/metas/base.hpp
  isVerificationFile: false
  path: cp-algo/data_structures/segment_tree/metas/affine.hpp
  requiredBy: []
  timestamp: '2024-02-10 22:49:03+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
documentation_of: cp-algo/data_structures/segment_tree/metas/affine.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/segment_tree/metas/affine.hpp
- /library/cp-algo/data_structures/segment_tree/metas/affine.hpp.html
title: cp-algo/data_structures/segment_tree/metas/affine.hpp
---
