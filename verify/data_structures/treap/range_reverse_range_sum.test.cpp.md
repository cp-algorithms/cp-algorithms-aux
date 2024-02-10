---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/affine.hpp
    title: cp-algo/algebra/affine.hpp
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
    path: cp-algo/data_structures/treap/metas/reverse.hpp
    title: cp-algo/data_structures/treap/metas/reverse.hpp
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
    PROBLEM: https://judge.yosupo.jp/problem/range_reverse_range_sum
    document_title: Range Reverse Range Sum
    links:
    - https://judge.yosupo.jp/problem/range_reverse_range_sum
  bundledCode: "#line 1 \"verify/data_structures/treap/range_reverse_range_sum.test.cpp\"\
    \n// @brief Range Reverse Range Sum\n#define PROBLEM \"https://judge.yosupo.jp/problem/range_reverse_range_sum\"\
    \n#line 1 \"cp-algo/data_structures/treap/metas/reverse.hpp\"\n\n\n#line 1 \"\
    cp-algo/data_structures/treap/metas/base.hpp\"\n\n\n#line 1 \"cp-algo/data_structures/treap/common.hpp\"\
    \n\n\n#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())\n\
    \n#line 4 \"cp-algo/data_structures/treap/metas/base.hpp\"\n#include <functional>\n\
    #include <algorithm>\n#include <cstdint>\n#define _safe_meta(i, op) _safe(i, _meta.op)\n\
    namespace cp_algo::data_structures::treap::metas {\n    struct base_meta {\n \
    \       void pull(auto const, auto const){}\n        void push(auto&, auto&){}\n\
    \    };\n}\n\n#line 1 \"cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n\
    #include <cassert>\nnamespace cp_algo::algebra {\n    template<typename base>\n\
    \    // a * x + b\n    struct lin {\n        base a = 1, b = 0;\n        std::optional<base>\
    \ c;\n        lin() {}\n        lin(base b): a(0), b(b) {}\n        lin(base a,\
    \ base b): a(a), b(b) {}\n        lin(base a, base b, base _c): a(a), b(b), c(_c)\
    \ {}\n\n        // polynomial product modulo x^2 - c\n        lin operator * (const\
    \ lin& t) {\n            assert(c && t.c && *c == *t.c);\n            return {a\
    \ * t.b + b * t.a, b * t.b + a * t.a * (*c), *c};\n        }\n\n        // a *\
    \ (t.a * x + t.b) + b\n        lin apply(lin const& t) const {\n            return\
    \ {a * t.a, a * t.b + b};\n        }\n\n        void prepend(lin const& t) {\n\
    \            *this = t.apply(*this);\n        }\n\n        base eval(base x) const\
    \ {\n            return a * x + b;\n        }\n    };\n}\n\n#line 6 \"cp-algo/data_structures/treap/metas/reverse.hpp\"\
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
    \                }\n            }\n        };\n}\n\n#line 1 \"cp-algo/data_structures/treap.hpp\"\
    \n\n\n#line 1 \"cp-algo/random/rng.hpp\"\n\n\n#include <chrono>\n#include <random>\n\
    namespace cp_algo::random {\n    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());\
    \ \n}\n\n#line 5 \"cp-algo/data_structures/treap.hpp\"\n#include <array>\n/* Submissions\
    \ on Library Judge:\n  Range Reverse Range Sum, 558ms - https://judge.yosupo.jp/submission/147860\n\
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
    }\n\n#line 5 \"verify/data_structures/treap/range_reverse_range_sum.test.cpp\"\
    \n#include <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::data_structures::treap;\n\
    \nusing meta = metas::reverse_meta<int64_t>;\nusing node = treap_node<meta>;\n\
    using treap = node::treap;\n\nvoid solve() {\n    istream_iterator<int> input(cin);\n\
    \    int n = *input++;\n    int q = *input++;\n    vector<treap> nodes(n);\n \
    \   generate_n(begin(nodes), n, [&](){\n        return new node(meta(*input++));\n\
    \    });\n    treap me = node::build(nodes);\n    while(q--) {\n        int t\
    \ = *input++;\n        int l = *input++;\n        int r = *input++;\n        if(t\
    \ == 0) {\n            node::exec_on_segment(me, l, r, [](treap &t) {\n      \
    \          _safe_meta(t, reverse = true);\n            });\n        } else {\n\
    \            node::exec_on_segment(me, l, r, [](treap const& t) {\n          \
    \      cout << _safe_meta(t, sum) << \"\\n\";\n            });\n        }\n  \
    \  }\n}\n\nsigned main() {\n    //freopen(\"input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n\
    \    cin.tie(0);\n    int t = 1;\n    while(t--) {\n        solve();\n    }\n\
    }\n"
  code: "// @brief Range Reverse Range Sum\n#define PROBLEM \"https://judge.yosupo.jp/problem/range_reverse_range_sum\"\
    \n#include \"cp-algo/data_structures/treap/metas/reverse.hpp\"\n#include \"cp-algo/data_structures/treap.hpp\"\
    \n#include <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::data_structures::treap;\n\
    \nusing meta = metas::reverse_meta<int64_t>;\nusing node = treap_node<meta>;\n\
    using treap = node::treap;\n\nvoid solve() {\n    istream_iterator<int> input(cin);\n\
    \    int n = *input++;\n    int q = *input++;\n    vector<treap> nodes(n);\n \
    \   generate_n(begin(nodes), n, [&](){\n        return new node(meta(*input++));\n\
    \    });\n    treap me = node::build(nodes);\n    while(q--) {\n        int t\
    \ = *input++;\n        int l = *input++;\n        int r = *input++;\n        if(t\
    \ == 0) {\n            node::exec_on_segment(me, l, r, [](treap &t) {\n      \
    \          _safe_meta(t, reverse = true);\n            });\n        } else {\n\
    \            node::exec_on_segment(me, l, r, [](treap const& t) {\n          \
    \      cout << _safe_meta(t, sum) << \"\\n\";\n            });\n        }\n  \
    \  }\n}\n\nsigned main() {\n    //freopen(\"input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n\
    \    cin.tie(0);\n    int t = 1;\n    while(t--) {\n        solve();\n    }\n\
    }\n"
  dependsOn:
  - cp-algo/data_structures/treap/metas/reverse.hpp
  - cp-algo/data_structures/treap/metas/base.hpp
  - cp-algo/data_structures/treap/common.hpp
  - cp-algo/algebra/affine.hpp
  - cp-algo/data_structures/treap.hpp
  - cp-algo/random/rng.hpp
  - cp-algo/data_structures/treap/common.hpp
  isVerificationFile: true
  path: verify/data_structures/treap/range_reverse_range_sum.test.cpp
  requiredBy: []
  timestamp: '2024-02-11 00:07:44+01:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/data_structures/treap/range_reverse_range_sum.test.cpp
layout: document
redirect_from:
- /verify/verify/data_structures/treap/range_reverse_range_sum.test.cpp
- /verify/verify/data_structures/treap/range_reverse_range_sum.test.cpp.html
title: Range Reverse Range Sum
---
