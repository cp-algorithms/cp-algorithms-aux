---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap.hpp
    title: cp-algo/data_structures/treap.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/common.hpp
    title: cp-algo/data_structures/treap/common.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/common.hpp
    title: cp-algo/data_structures/treap/common.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/treap/metas/base.hpp
    title: cp-algo/data_structures/treap/metas/base.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/random/rng.hpp
    title: cp-algo/random/rng.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/cartesian_tree
    document_title: Build Cartesian Tree
    links:
    - https://judge.yosupo.jp/problem/cartesian_tree
  bundledCode: "#line 1 \"verify/data_structures/treap/cartesian_tree.test.cpp\"\n\
    // @brief Build Cartesian Tree\n#define PROBLEM \"https://judge.yosupo.jp/problem/cartesian_tree\"\
    \n#line 1 \"cp-algo/data_structures/treap/metas/base.hpp\"\n\n\n#line 1 \"cp-algo/data_structures/treap/common.hpp\"\
    \n\n\n#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())\n\
    \n#line 4 \"cp-algo/data_structures/treap/metas/base.hpp\"\n#include <functional>\n\
    #include <algorithm>\n#include <cstdint>\n#define _safe_meta(i, op) _safe(i, _meta.op)\n\
    namespace cp_algo::data_structures::treap::metas {\n    struct base_meta {\n \
    \       void pull(auto const, auto const){}\n        void push(auto&, auto&){}\n\
    \    };\n}\n\n#line 1 \"cp-algo/data_structures/treap.hpp\"\n\n\n#line 1 \"cp-algo/random/rng.hpp\"\
    \n\n\n#include <chrono>\n#include <random>\nnamespace cp_algo::random {\n    std::mt19937_64\
    \ rng(std::chrono::steady_clock::now().time_since_epoch().count()); \n}\n\n#line\
    \ 5 \"cp-algo/data_structures/treap.hpp\"\n#include <array>\n/* Submissions on\
    \ Library Judge:\n  Range Reverse Range Sum, 558ms - https://judge.yosupo.jp/submission/147860\n\
    \  Cartesian Tree, 229ms - https://judge.yosupo.jp/submission/147858\n  Dynamic\
    \ Sequence Range Affine Range Sum, 2245ms - https://judge.yosupo.jp/submission/148948\n\
    */\nnamespace cp_algo::data_structures::treap {\n    template<typename meta>\n\
    \    struct treap_node {\n\n        using node = treap_node;\n        using treap\
    \ = node*;\n        meta _meta;\n        int prior = random::rng();\n        size_t\
    \ size = 1;\n        treap children[2] = {nullptr, nullptr};\n        enum subtree\
    \ {L, R};\n\n        treap pull() {\n            _meta.pull(children[L], children[R]);\n\
    \            size = 1 + _safe(children[L], size) + _safe(children[R], size);\n\
    \            return this;\n        }\n\n        treap push() {\n            _meta.push(children[L],\
    \ children[R]);\n            return this;\n        }\n\n        // set i-th child\
    \ and pull metadata\n        treap set(subtree i, treap t) {\n            children[i]\
    \ = t;\n            return pull();\n        }\n\n        // push changes and detach\
    \ the i-th child\n        treap cut(subtree i) {\n            return children[i];\n\
    \        }\n\n        static treap merge(treap A, treap B) {\n            if(!_safe(A,\
    \ push()) || !_safe(B, push())) {\n                return A ? A : B;\n       \
    \     } else if(A->prior < B->prior) {\n                return A->set(R, merge(A->cut(R),\
    \ B));\n            } else {\n                return B->set(L, merge(A, B->cut(L)));\n\
    \            }\n        }\n\n        // return {L, R}, where |L|=k or L=A when\
    \ |A| < k\n        static std::array<treap, 2> split(treap A, size_t k) {\n  \
    \          if(!_safe(A, push())) {\n                return {nullptr, nullptr};\n\
    \            } else if(_safe(A->children[L], size) >= k) {\n                auto\
    \ [split_L, split_R] = split(A->cut(L), k);\n                return {split_L,\
    \ A->set(L, split_R)};\n            } else {\n                k -= _safe(A->children[L],\
    \ size) + 1;\n                auto [split_L, split_R] = split(A->cut(R), k);\n\
    \                return {A->set(R, split_L), split_R};\n            }\n      \
    \  }\n\n        static void exec_on_segment(treap &A, size_t l, size_t r, auto\
    \ func) {\n            auto [LM, R] = split(A, r);\n            auto [L, M] =\
    \ split(LM, l);\n            func(M);\n            A = merge(L, merge(M, R));\n\
    \        }\n\n        static void insert(treap &A, size_t pos, treap t) {\n  \
    \          auto [L, R] = split(A, pos);\n            A = merge(L, merge(t, R));\n\
    \        }\n\n        static void erase(treap &A, size_t pos) {\n            auto\
    \ [L, MR] = split(A, pos);\n            auto [M, R] = split(MR, 1);\n        \
    \    delete M;\n            A = merge(L, R);\n        }\n\n        static void\
    \ exec_on_each(treap &A, auto func) {\n            if(A) {\n                exec_on_each(A->children[L],\
    \ func);\n                func(A);\n                exec_on_each(A->children[R],\
    \ func);\n            }\n        }\n\n        treap pull_all() {\n           \
    \ _safe(children[L], pull_all());\n            _safe(children[R], pull_all());\n\
    \            return pull();\n        }\n\n        treap push_all() {\n       \
    \     push();\n            _safe(children[L], push_all());\n            _safe(children[R],\
    \ push_all());\n            return this;\n        }\n\n        static treap build(auto\
    \ const& nodes) {\n            std::vector<treap> st;\n            for(auto cur:\
    \ nodes) {\n                while(st.size() >= 2 && st[st.size() - 2]->prior >\
    \ cur->prior) {\n                    st.pop_back();\n                }\n     \
    \           if(!st.empty() && st.back()->prior > cur->prior) {\n             \
    \       cur->set(L, st.back());\n                    st.pop_back();\n        \
    \        }\n                if(!st.empty() && st.back()->prior < cur->prior) {\n\
    \                    st.back()->set(R, cur);\n                }\n            \
    \    st.push_back(cur);\n            }\n            return st.empty() ? nullptr\
    \ : st[0]->pull_all();\n        }\n    };\n\n    struct null_meta {\n        void\
    \ pull(auto const, auto const) {}\n        void push(auto&, auto&) {}\n    };\n\
    }\n\n#line 5 \"verify/data_structures/treap/cartesian_tree.test.cpp\"\n#include\
    \ <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::data_structures::treap;\n\
    \nstruct val_meta: metas::base_meta {\n    int val;\n    val_meta(int val): val(val){}\n\
    };\n\nusing node = treap_node<val_meta>;\nusing treap = node::treap;\n\nvoid solve()\
    \ {\n    istream_iterator<int> input(cin);\n    int n = *input++;\n    vector<treap>\
    \ nodes(n);\n    for(int i = 0; i < n; i++) {\n        nodes[i] = new node(val_meta(i),\
    \ *input++);\n    }\n    treap me = node::build(nodes);\n    vector<int> p(n,\
    \ -1);\n    node::exec_on_each(me, [&](treap t) {\n        for(auto child: t->children)\
    \ {\n            if(child) {\n                p[_safe_meta(child, val)] = _safe_meta(t,\
    \ val);\n            }\n        }\n    });\n    for(int i = 0; i < n; i++) {\n\
    \        cout << (p[i] == -1 ? i : p[i]) << ' ';\n    }\n}\n\nsigned main() {\n\
    \    //freopen(\"input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n \
    \   cin.tie(0);\n    int t = 1;\n    while(t--) {\n        solve();\n    }\n}\n"
  code: "// @brief Build Cartesian Tree\n#define PROBLEM \"https://judge.yosupo.jp/problem/cartesian_tree\"\
    \n#include \"cp-algo/data_structures/treap/metas/base.hpp\"\n#include \"cp-algo/data_structures/treap.hpp\"\
    \n#include <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::data_structures::treap;\n\
    \nstruct val_meta: metas::base_meta {\n    int val;\n    val_meta(int val): val(val){}\n\
    };\n\nusing node = treap_node<val_meta>;\nusing treap = node::treap;\n\nvoid solve()\
    \ {\n    istream_iterator<int> input(cin);\n    int n = *input++;\n    vector<treap>\
    \ nodes(n);\n    for(int i = 0; i < n; i++) {\n        nodes[i] = new node(val_meta(i),\
    \ *input++);\n    }\n    treap me = node::build(nodes);\n    vector<int> p(n,\
    \ -1);\n    node::exec_on_each(me, [&](treap t) {\n        for(auto child: t->children)\
    \ {\n            if(child) {\n                p[_safe_meta(child, val)] = _safe_meta(t,\
    \ val);\n            }\n        }\n    });\n    for(int i = 0; i < n; i++) {\n\
    \        cout << (p[i] == -1 ? i : p[i]) << ' ';\n    }\n}\n\nsigned main() {\n\
    \    //freopen(\"input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n \
    \   cin.tie(0);\n    int t = 1;\n    while(t--) {\n        solve();\n    }\n}\n"
  dependsOn:
  - cp-algo/data_structures/treap/metas/base.hpp
  - cp-algo/data_structures/treap/common.hpp
  - cp-algo/data_structures/treap.hpp
  - cp-algo/random/rng.hpp
  - cp-algo/data_structures/treap/common.hpp
  isVerificationFile: true
  path: verify/data_structures/treap/cartesian_tree.test.cpp
  requiredBy: []
  timestamp: '2024-02-11 11:53:49+01:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/data_structures/treap/cartesian_tree.test.cpp
layout: document
redirect_from:
- /verify/verify/data_structures/treap/cartesian_tree.test.cpp
- /verify/verify/data_structures/treap/cartesian_tree.test.cpp.html
title: Build Cartesian Tree
---
