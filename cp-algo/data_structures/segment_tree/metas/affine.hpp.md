---
data:
  _extendedDependsOn:
  - icon: ':question:'
    path: cp-algo/algebra/affine.hpp
    title: cp-algo/algebra/affine.hpp
  - icon: ':x:'
    path: cp-algo/data_structures/segment_tree/metas/base.hpp
    title: cp-algo/data_structures/segment_tree/metas/base.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':x:'
    path: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
    title: Range Affine Range Sum
  _isVerificationFailed: true
  _pathExtension: hpp
  _verificationStatusIcon: ':x:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/data_structures/segment_tree/metas/affine.hpp\"\n\
    \n\n#line 1 \"cp-algo/data_structures/segment_tree/metas/base.hpp\"\n\n\n#include\
    \ <functional>\n#include <algorithm>\n#include <cstdint>\nnamespace cp_algo::data_structures::segment_tree::metas\
    \ {\n    template<typename derived_meta>\n    struct base_meta {\n        using\
    \ meta = derived_meta;\n        void pull(meta const&, meta const&, int, int)\
    \ {};\n        void push(meta*, meta*, int, int) {};\n    };\n}\n\n#line 1 \"\
    cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n#include <cassert>\nnamespace\
    \ cp_algo::algebra {\n    template<typename base>\n    // a * x + b\n    struct\
    \ lin {\n        base a = 1, b = 0;\n        std::optional<base> c;\n        lin()\
    \ {}\n        lin(base b): a(0), b(b) {}\n        lin(base a, base b): a(a), b(b)\
    \ {}\n        lin(base a, base b, base _c): a(a), b(b), c(_c) {}\n\n        //\
    \ polynomial product modulo x^2 - c\n        lin operator * (const lin& t) {\n\
    \            assert(c && t.c && *c == *t.c);\n            return lin(a * t.b +\
    \ b * t.a, b * t.b + a * t.a * (*c), *c);\n        }\n\n        // a * (t.a *\
    \ x + t.b) + b\n        lin compose(lin const& t) const {\n            return\
    \ lin{a * t.a, a * t.b + b};\n        }\n\n        void prepend(lin const& t)\
    \ {\n            *this = t.compose(*this);\n        }\n\n        base eval(base\
    \ x) const {\n            return a * x + b;\n        }\n    };\n}\n\n#line 5 \"\
    cp-algo/data_structures/segment_tree/metas/affine.hpp\"\nnamespace cp_algo::data_structures::segment_tree::metas\
    \ {\n    template<typename base>\n    struct affine_meta: base_meta<affine_meta<base>>\
    \ {\n        using meta = affine_meta;\n\n        base sum = 0;\n        algebra<base>::lin\
    \ to_push = {};\n\n        affine_meta() {}\n        affine_meta(base sum): sum(sum)\
    \ {}\n\n        void push(meta *L, meta *R, int l, int r) override {\n       \
    \     if(to_push.a != 1 || to_push.b != 0) {\n                sum = to_push.a\
    \ * sum + to_push.b * (r - l);\n                if(r - l > 1) {\n            \
    \        L->to_push = to_push.compose(L->to_push);\n                    R->to_push\
    \ = to_push.compose(R->to_push);\n                }\n                to_push =\
    \ {};\n            }\n        }\n\n        void pull(meta const& L, meta const&\
    \ R, int, int) override {\n            sum = L.sum + R.sum;\n        }\n    };\n\
    }\n\n"
  code: "#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP\n#define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP\n\
    #include \"base.hpp\"\n#include \"../../../algebra/affine.hpp\"\nnamespace cp_algo::data_structures::segment_tree::metas\
    \ {\n    template<typename base>\n    struct affine_meta: base_meta<affine_meta<base>>\
    \ {\n        using meta = affine_meta;\n\n        base sum = 0;\n        algebra<base>::lin\
    \ to_push = {};\n\n        affine_meta() {}\n        affine_meta(base sum): sum(sum)\
    \ {}\n\n        void push(meta *L, meta *R, int l, int r) override {\n       \
    \     if(to_push.a != 1 || to_push.b != 0) {\n                sum = to_push.a\
    \ * sum + to_push.b * (r - l);\n                if(r - l > 1) {\n            \
    \        L->to_push = to_push.compose(L->to_push);\n                    R->to_push\
    \ = to_push.compose(R->to_push);\n                }\n                to_push =\
    \ {};\n            }\n        }\n\n        void pull(meta const& L, meta const&\
    \ R, int, int) override {\n            sum = L.sum + R.sum;\n        }\n    };\n\
    }\n#endif // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP"
  dependsOn:
  - cp-algo/data_structures/segment_tree/metas/base.hpp
  - cp-algo/algebra/affine.hpp
  isVerificationFile: false
  path: cp-algo/data_structures/segment_tree/metas/affine.hpp
  requiredBy: []
  timestamp: '2024-02-10 23:55:00+01:00'
  verificationStatus: LIBRARY_ALL_WA
  verifiedWith:
  - verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
documentation_of: cp-algo/data_structures/segment_tree/metas/affine.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/segment_tree/metas/affine.hpp
- /library/cp-algo/data_structures/segment_tree/metas/affine.hpp.html
title: cp-algo/data_structures/segment_tree/metas/affine.hpp
---
