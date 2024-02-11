---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/affine.hpp
    title: cp-algo/algebra/affine.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/common.hpp
    title: cp-algo/data_structures/treap/common.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/metas/base.hpp
    title: cp-algo/data_structures/treap/metas/base.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith:
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
  bundledCode: "#line 1 \"cp-algo/data_structures/treap/metas/reverse.hpp\"\n\n\n\
    #line 1 \"cp-algo/data_structures/treap/metas/base.hpp\"\n\n\n#line 1 \"cp-algo/data_structures/treap/common.hpp\"\
    \n\n\n#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())\n\
    \n#line 4 \"cp-algo/data_structures/treap/metas/base.hpp\"\n#include <functional>\n\
    #include <algorithm>\n#include <cstdint>\n#define _safe_meta(i, op) _safe(i, _meta.op)\n\
    namespace cp_algo::data_structures::treap::metas {\n    struct base_meta {\n \
    \       void pull(auto const, auto const){}\n        void push(auto&, auto&){}\n\
    \    };\n}\n\n#line 1 \"cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n\
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
    \ * A + b * B, c * A + d * B};\n        }\n    };\n}\n\n#line 6 \"cp-algo/data_structures/treap/metas/reverse.hpp\"\
    \nnamespace cp_algo::data_structures::treap::metas {\n        template<typename\
    \ base>\n        struct reverse_meta: base_meta {\n            using lin = algebra::lin<base>;\n\
    \            base val;\n            size_t sz = 1;\n            bool reverse =\
    \ false;\n            base sum = val;\n            \n            lin to_push =\
    \ {};\n\n            reverse_meta(base val): val(val) {}\n\n            void pull(auto\
    \ const L, auto const R) {\n                sum = val + _safe_meta(L, sum) + _safe_meta(R,\
    \ sum);\n                sz = 1 + _safe_meta(L, sz) + _safe_meta(R, sz);\n   \
    \         }\n            void add_push(lin const& t) {\n                val =\
    \ t.eval(val);\n                sum = t.a * sum + t.b * sz;\n                to_push.prepend(t);\n\
    \            }\n            void push(auto &L, auto &R) {\n                if(reverse)\
    \ {\n                    reverse = false;\n                    std::swap(L, R);\n\
    \                    _safe_meta(L, reverse ^= 1);\n                    _safe_meta(R,\
    \ reverse ^= 1);\n                }\n                if(to_push.a != 1 || to_push.b\
    \ != 0) {\n                    _safe_meta(L, add_push(to_push));\n           \
    \         _safe_meta(R, add_push(to_push));\n                    to_push = {};\n\
    \                }\n            }\n        };\n}\n\n"
  code: "#ifndef CP_ALGO_DATA_STRUCTURES_TREAP_METAS_REVERSE_HPP\n#define CP_ALGO_DATA_STRUCTURES_TREAP_METAS_REVERSE_HPP\n\
    #include \"base.hpp\"\n#include \"../../../algebra/affine.hpp\"\n#include <algorithm>\n\
    namespace cp_algo::data_structures::treap::metas {\n        template<typename\
    \ base>\n        struct reverse_meta: base_meta {\n            using lin = algebra::lin<base>;\n\
    \            base val;\n            size_t sz = 1;\n            bool reverse =\
    \ false;\n            base sum = val;\n            \n            lin to_push =\
    \ {};\n\n            reverse_meta(base val): val(val) {}\n\n            void pull(auto\
    \ const L, auto const R) {\n                sum = val + _safe_meta(L, sum) + _safe_meta(R,\
    \ sum);\n                sz = 1 + _safe_meta(L, sz) + _safe_meta(R, sz);\n   \
    \         }\n            void add_push(lin const& t) {\n                val =\
    \ t.eval(val);\n                sum = t.a * sum + t.b * sz;\n                to_push.prepend(t);\n\
    \            }\n            void push(auto &L, auto &R) {\n                if(reverse)\
    \ {\n                    reverse = false;\n                    std::swap(L, R);\n\
    \                    _safe_meta(L, reverse ^= 1);\n                    _safe_meta(R,\
    \ reverse ^= 1);\n                }\n                if(to_push.a != 1 || to_push.b\
    \ != 0) {\n                    _safe_meta(L, add_push(to_push));\n           \
    \         _safe_meta(R, add_push(to_push));\n                    to_push = {};\n\
    \                }\n            }\n        };\n}\n#endif // CP_ALGO_DATA_STRUCTURES_TREAP_METAS_REVERSE_HPP"
  dependsOn:
  - cp-algo/data_structures/treap/metas/base.hpp
  - cp-algo/data_structures/treap/common.hpp
  - cp-algo/algebra/affine.hpp
  isVerificationFile: false
  path: cp-algo/data_structures/treap/metas/reverse.hpp
  requiredBy: []
  timestamp: '2024-02-11 13:25:01+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/treap/range_reverse_range_sum.test.cpp
  - verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
documentation_of: cp-algo/data_structures/treap/metas/reverse.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/treap/metas/reverse.hpp
- /library/cp-algo/data_structures/treap/metas/reverse.hpp.html
title: cp-algo/data_structures/treap/metas/reverse.hpp
---
