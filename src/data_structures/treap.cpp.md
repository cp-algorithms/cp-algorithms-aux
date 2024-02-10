---
data:
  _extendedDependsOn: []
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':warning:'
  attributes:
    links:
    - https://judge.yosupo.jp/submission/147858
    - https://judge.yosupo.jp/submission/147860
    - https://judge.yosupo.jp/submission/148948
  bundledCode: "#line 1 \"src/data_structures/treap.cpp\"\n/* Submissions on Library\
    \ Judge:\n  Range Reverse Range Sum, 558ms - https://judge.yosupo.jp/submission/147860\n\
    \  Cartesian Tree, 229ms - https://judge.yosupo.jp/submission/147858\n  Dynamic\
    \ Sequence Range Affine Range Sum, 2245ms - https://judge.yosupo.jp/submission/148948\n\
    \  */\nnamespace data_structures {\n    namespace treap {\n        mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());\n\
    \        #define safe(t, op) (t ? t->op : typename remove_reference<decltype(t->op)>::type())\n\
    \        template<typename meta>\n        struct treap_node {\n\n            using\
    \ node = treap_node;\n            using treap = node*;\n            meta _meta;\n\
    \            int prior = rng();\n            size_t size = 1;\n            treap\
    \ children[2] = {nullptr, nullptr};\n            enum subtree {L, R};\n\n    \
    \        base check_sum() {\n                push();\n                return _meta.val\
    \ + safe(children[L], check_sum()) + safe(children[R], check_sum());\n       \
    \     }\n\n            treap pull() {\n                _meta.pull(children[L],\
    \ children[R]);\n                size = 1 + safe(children[L], size) + safe(children[R],\
    \ size);\n                return this;\n            }\n\n            treap push()\
    \ {\n                _meta.push(children[L], children[R]);\n                return\
    \ this;\n            }\n\n            // set i-th child and pull metadata\n  \
    \          treap set(subtree i, treap t) {\n                children[i] = t;\n\
    \                return pull();\n            }\n\n            // push changes\
    \ and detach the i-th child\n            treap cut(subtree i) {\n            \
    \    return children[i];\n            }\n\n            static treap merge(treap\
    \ A, treap B) {\n                if(!safe(A, push()) || !safe(B, push())) {\n\
    \                    return A ? A : B;\n                } else if(A->prior < B->prior)\
    \ {\n                    return A->set(R, merge(A->cut(R), B));\n            \
    \    } else {\n                    return B->set(L, merge(A, B->cut(L)));\n  \
    \              }\n            }\n\n            // return {L, R}, where |L|=k or\
    \ L=A when |A| < k\n            static pair<treap, treap> split(treap A, size_t\
    \ k) {\n                if(!safe(A, push())) {\n                    return {nullptr,\
    \ nullptr};\n                } else if(safe(A->children[L], size) >= k) {\n  \
    \                  auto [split_L, split_R] = split(A->cut(L), k);\n          \
    \          return {split_L, A->set(L, split_R)};\n                } else {\n \
    \                   k -= safe(A->children[L], size) + 1;\n                   \
    \ auto [split_L, split_R] = split(A->cut(R), k);\n                    return {A->set(R,\
    \ split_L), split_R};\n                }\n            }\n\n            static\
    \ void exec_on_segment(treap &A, size_t l, size_t r, auto func) {\n          \
    \      auto [LM, R] = split(A, r);\n                auto [L, M] = split(LM, l);\n\
    \                func(M);\n                A = merge(L, merge(M, R));\n      \
    \      }\n\n            static void insert(treap &A, size_t pos, treap t) {\n\
    \                auto [L, R] = split(A, pos);\n                A = merge(L, merge(t,\
    \ R));\n            }\n\n            static void erase(treap &A, size_t pos) {\n\
    \                auto [L, MR] = split(A, pos);\n                auto [M, R] =\
    \ split(MR, 1);\n                delete M;\n                A = merge(L, R);\n\
    \            }\n\n            static void exec_on_each(treap &A, auto func) {\n\
    \                if(A) {\n                    exec_on_each(A->children[L], func);\n\
    \                    func(A);\n                    exec_on_each(A->children[R],\
    \ func);\n                }\n            }\n\n            treap pull_all() {\n\
    \                safe(children[L], pull_all());\n                safe(children[R],\
    \ pull_all());\n                return pull();\n            }\n\n            treap\
    \ push_all() {\n                push();\n                safe(children[L], push_all());\n\
    \                safe(children[R], push_all());\n                return this;\n\
    \            }\n\n            static treap build(vector<treap> nodes) {\n    \
    \            vector<treap> st;\n                for(auto cur: nodes) {\n     \
    \               while(st.size() >= 2 && st[st.size() - 2]->prior > cur->prior)\
    \ {\n                        st.pop_back();\n                    }\n         \
    \           if(!st.empty() && st.back()->prior > cur->prior) {\n             \
    \           cur->set(L, st.back());\n                        st.pop_back();\n\
    \                    }\n                    if(!st.empty() && st.back()->prior\
    \ < cur->prior) {\n                        st.back()->set(R, cur);\n         \
    \           }\n                    st.push_back(cur);\n                }\n   \
    \             return st.empty() ? nullptr : st[0]->pull_all();\n            }\n\
    \        };\n\n        struct null_meta {\n            void pull(auto const, auto\
    \ const) {}\n            void push(auto&, auto&) {}\n        };\n    }\n}\n"
  code: "/* Submissions on Library Judge:\n  Range Reverse Range Sum, 558ms - https://judge.yosupo.jp/submission/147860\n\
    \  Cartesian Tree, 229ms - https://judge.yosupo.jp/submission/147858\n  Dynamic\
    \ Sequence Range Affine Range Sum, 2245ms - https://judge.yosupo.jp/submission/148948\n\
    \  */\nnamespace data_structures {\n    namespace treap {\n        mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());\n\
    \        #define safe(t, op) (t ? t->op : typename remove_reference<decltype(t->op)>::type())\n\
    \        template<typename meta>\n        struct treap_node {\n\n            using\
    \ node = treap_node;\n            using treap = node*;\n            meta _meta;\n\
    \            int prior = rng();\n            size_t size = 1;\n            treap\
    \ children[2] = {nullptr, nullptr};\n            enum subtree {L, R};\n\n    \
    \        base check_sum() {\n                push();\n                return _meta.val\
    \ + safe(children[L], check_sum()) + safe(children[R], check_sum());\n       \
    \     }\n\n            treap pull() {\n                _meta.pull(children[L],\
    \ children[R]);\n                size = 1 + safe(children[L], size) + safe(children[R],\
    \ size);\n                return this;\n            }\n\n            treap push()\
    \ {\n                _meta.push(children[L], children[R]);\n                return\
    \ this;\n            }\n\n            // set i-th child and pull metadata\n  \
    \          treap set(subtree i, treap t) {\n                children[i] = t;\n\
    \                return pull();\n            }\n\n            // push changes\
    \ and detach the i-th child\n            treap cut(subtree i) {\n            \
    \    return children[i];\n            }\n\n            static treap merge(treap\
    \ A, treap B) {\n                if(!safe(A, push()) || !safe(B, push())) {\n\
    \                    return A ? A : B;\n                } else if(A->prior < B->prior)\
    \ {\n                    return A->set(R, merge(A->cut(R), B));\n            \
    \    } else {\n                    return B->set(L, merge(A, B->cut(L)));\n  \
    \              }\n            }\n\n            // return {L, R}, where |L|=k or\
    \ L=A when |A| < k\n            static pair<treap, treap> split(treap A, size_t\
    \ k) {\n                if(!safe(A, push())) {\n                    return {nullptr,\
    \ nullptr};\n                } else if(safe(A->children[L], size) >= k) {\n  \
    \                  auto [split_L, split_R] = split(A->cut(L), k);\n          \
    \          return {split_L, A->set(L, split_R)};\n                } else {\n \
    \                   k -= safe(A->children[L], size) + 1;\n                   \
    \ auto [split_L, split_R] = split(A->cut(R), k);\n                    return {A->set(R,\
    \ split_L), split_R};\n                }\n            }\n\n            static\
    \ void exec_on_segment(treap &A, size_t l, size_t r, auto func) {\n          \
    \      auto [LM, R] = split(A, r);\n                auto [L, M] = split(LM, l);\n\
    \                func(M);\n                A = merge(L, merge(M, R));\n      \
    \      }\n\n            static void insert(treap &A, size_t pos, treap t) {\n\
    \                auto [L, R] = split(A, pos);\n                A = merge(L, merge(t,\
    \ R));\n            }\n\n            static void erase(treap &A, size_t pos) {\n\
    \                auto [L, MR] = split(A, pos);\n                auto [M, R] =\
    \ split(MR, 1);\n                delete M;\n                A = merge(L, R);\n\
    \            }\n\n            static void exec_on_each(treap &A, auto func) {\n\
    \                if(A) {\n                    exec_on_each(A->children[L], func);\n\
    \                    func(A);\n                    exec_on_each(A->children[R],\
    \ func);\n                }\n            }\n\n            treap pull_all() {\n\
    \                safe(children[L], pull_all());\n                safe(children[R],\
    \ pull_all());\n                return pull();\n            }\n\n            treap\
    \ push_all() {\n                push();\n                safe(children[L], push_all());\n\
    \                safe(children[R], push_all());\n                return this;\n\
    \            }\n\n            static treap build(vector<treap> nodes) {\n    \
    \            vector<treap> st;\n                for(auto cur: nodes) {\n     \
    \               while(st.size() >= 2 && st[st.size() - 2]->prior > cur->prior)\
    \ {\n                        st.pop_back();\n                    }\n         \
    \           if(!st.empty() && st.back()->prior > cur->prior) {\n             \
    \           cur->set(L, st.back());\n                        st.pop_back();\n\
    \                    }\n                    if(!st.empty() && st.back()->prior\
    \ < cur->prior) {\n                        st.back()->set(R, cur);\n         \
    \           }\n                    st.push_back(cur);\n                }\n   \
    \             return st.empty() ? nullptr : st[0]->pull_all();\n            }\n\
    \        };\n\n        struct null_meta {\n            void pull(auto const, auto\
    \ const) {}\n            void push(auto&, auto&) {}\n        };\n    }\n}\n"
  dependsOn: []
  isVerificationFile: false
  path: src/data_structures/treap.cpp
  requiredBy: []
  timestamp: '2024-02-10 16:40:11+01:00'
  verificationStatus: LIBRARY_NO_TESTS
  verifiedWith: []
documentation_of: src/data_structures/treap.cpp
layout: document
redirect_from:
- /library/src/data_structures/treap.cpp
- /library/src/data_structures/treap.cpp.html
title: src/data_structures/treap.cpp
---
