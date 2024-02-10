---
data:
  _extendedDependsOn: []
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':warning:'
  attributes:
    links:
    - https://judge.yosupo.jp/submission/148324
    - https://judge.yosupo.jp/submission/148544
  bundledCode: "#line 1 \"cp-algo/data_structures/segment_tree.hpp\"\n\n\n#include\
    \ <vector>\n/* Metas examples:\n    Range Affine Range Sum, 1138ms - https://judge.yosupo.jp/submission/148324\n\
    \    Range Chmin Chmax Add Range Sum, 787ms - https://judge.yosupo.jp/submission/148544\n\
    */\nnamespace data_structures {\n    namespace segment_tree {\n        template<typename\
    \ meta>\n        struct segment_tree {\n            const int N;\n           \
    \ std::vector<meta> _meta;\n\n            segment_tree(int n): N(n), _meta(4 *\
    \ N) {}\n\n            segment_tree(std::vector<meta> leafs): N(size(leafs)),\
    \ _meta(4 * N) {\n                build(leafs);\n            }\n\n           \
    \ void pull(int v, int l, int r) {\n                if(r - l > 1) {\n        \
    \            _meta[v].pull(_meta[2 * v], _meta[2 * v + 1], l, r);\n          \
    \      }\n            }\n\n            void push(int v, int l, int r) {\n    \
    \            if(r - l > 1) {\n                    _meta[v].push(&_meta[2 * v],\
    \ &_meta[2 * v + 1], l, r);\n                } else {\n                    _meta[v].push(nullptr,\
    \ nullptr, l, r);\n                }\n            }\n\n            void build(auto\
    \ &a, int v, size_t l, size_t r) {\n                if(r - l == 1) {\n       \
    \             if(l < size(a)) {\n                        _meta[v] = a[l];\n  \
    \                  }\n                } else {\n                    size_t m =\
    \ (l + r) / 2;\n                    build(a, 2 * v, l, m);\n                 \
    \   build(a, 2 * v + 1, m, r);\n                    pull(v, l, r);\n         \
    \       }\n            }\n\n            void build(auto &a) {\n              \
    \  build(a, 1, 0, N);\n            }\n\n            void exec_on_segment(int a,\
    \ int b, auto func, auto proceed, auto stop, int v, int l, int r) {\n        \
    \        push(v, l, r);\n                if(r <= a || b <= l || stop(_meta[v]))\
    \ {\n                    return;\n                } else if(a <= l && r <= b &&\
    \ proceed(_meta[v])) {\n                    func(_meta[v]);\n                \
    \    push(v, l, r);\n                } else {\n                    int m = (l\
    \ + r) / 2;\n                    exec_on_segment(a, b, func, proceed, stop, 2\
    \ * v, l, m);\n                    exec_on_segment(a, b, func, proceed, stop,\
    \ 2 * v + 1, m, r);\n                    pull(v, l, r);\n                }\n \
    \           }\n\n            static constexpr auto default_true = [](auto const&){return\
    \ true;};\n            static constexpr auto default_false = [](auto const&){return\
    \ false;};\n\n            void exec_on_segment(int a, int b, auto func, auto proceed,\
    \ auto stop) {\n                exec_on_segment(a, b, func, proceed, stop, 1,\
    \ 0, N);\n            }\n\n            void exec_on_segment(int a, int b, auto\
    \ func) {\n                exec_on_segment(a, b, func, default_true, default_false);\n\
    \            }\n        };\n\n        namespace metas {\n            struct null_meta\
    \ {\n                using meta = null_meta;\n                void pull(meta const&,\
    \ meta const&, int, int) {}\n                void push(meta*, meta*, int, int)\
    \ {}\n            };\n        }\n    }\n}\n\n"
  code: "#ifndef DATA_STRUCTURES_SEGMENT_TREE_HPP\n#define DATA_STRUCTURES_SEGMENT_TREE_HPP\n\
    #include <vector>\n/* Metas examples:\n    Range Affine Range Sum, 1138ms - https://judge.yosupo.jp/submission/148324\n\
    \    Range Chmin Chmax Add Range Sum, 787ms - https://judge.yosupo.jp/submission/148544\n\
    */\nnamespace data_structures {\n    namespace segment_tree {\n        template<typename\
    \ meta>\n        struct segment_tree {\n            const int N;\n           \
    \ std::vector<meta> _meta;\n\n            segment_tree(int n): N(n), _meta(4 *\
    \ N) {}\n\n            segment_tree(std::vector<meta> leafs): N(size(leafs)),\
    \ _meta(4 * N) {\n                build(leafs);\n            }\n\n           \
    \ void pull(int v, int l, int r) {\n                if(r - l > 1) {\n        \
    \            _meta[v].pull(_meta[2 * v], _meta[2 * v + 1], l, r);\n          \
    \      }\n            }\n\n            void push(int v, int l, int r) {\n    \
    \            if(r - l > 1) {\n                    _meta[v].push(&_meta[2 * v],\
    \ &_meta[2 * v + 1], l, r);\n                } else {\n                    _meta[v].push(nullptr,\
    \ nullptr, l, r);\n                }\n            }\n\n            void build(auto\
    \ &a, int v, size_t l, size_t r) {\n                if(r - l == 1) {\n       \
    \             if(l < size(a)) {\n                        _meta[v] = a[l];\n  \
    \                  }\n                } else {\n                    size_t m =\
    \ (l + r) / 2;\n                    build(a, 2 * v, l, m);\n                 \
    \   build(a, 2 * v + 1, m, r);\n                    pull(v, l, r);\n         \
    \       }\n            }\n\n            void build(auto &a) {\n              \
    \  build(a, 1, 0, N);\n            }\n\n            void exec_on_segment(int a,\
    \ int b, auto func, auto proceed, auto stop, int v, int l, int r) {\n        \
    \        push(v, l, r);\n                if(r <= a || b <= l || stop(_meta[v]))\
    \ {\n                    return;\n                } else if(a <= l && r <= b &&\
    \ proceed(_meta[v])) {\n                    func(_meta[v]);\n                \
    \    push(v, l, r);\n                } else {\n                    int m = (l\
    \ + r) / 2;\n                    exec_on_segment(a, b, func, proceed, stop, 2\
    \ * v, l, m);\n                    exec_on_segment(a, b, func, proceed, stop,\
    \ 2 * v + 1, m, r);\n                    pull(v, l, r);\n                }\n \
    \           }\n\n            static constexpr auto default_true = [](auto const&){return\
    \ true;};\n            static constexpr auto default_false = [](auto const&){return\
    \ false;};\n\n            void exec_on_segment(int a, int b, auto func, auto proceed,\
    \ auto stop) {\n                exec_on_segment(a, b, func, proceed, stop, 1,\
    \ 0, N);\n            }\n\n            void exec_on_segment(int a, int b, auto\
    \ func) {\n                exec_on_segment(a, b, func, default_true, default_false);\n\
    \            }\n        };\n\n        namespace metas {\n            struct null_meta\
    \ {\n                using meta = null_meta;\n                void pull(meta const&,\
    \ meta const&, int, int) {}\n                void push(meta*, meta*, int, int)\
    \ {}\n            };\n        }\n    }\n}\n#endif // DATA_STRUCTURES_SEGMENT_TREE_HPP"
  dependsOn: []
  isVerificationFile: false
  path: cp-algo/data_structures/segment_tree.hpp
  requiredBy: []
  timestamp: '2024-02-10 19:59:38+01:00'
  verificationStatus: LIBRARY_NO_TESTS
  verifiedWith: []
documentation_of: cp-algo/data_structures/segment_tree.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/segment_tree.hpp
- /library/cp-algo/data_structures/segment_tree.hpp.html
title: cp-algo/data_structures/segment_tree.hpp
---
