---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/affine.hpp
    title: cp-algo/algebra/affine.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/affine.hpp
    title: cp-algo/algebra/affine.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/common.hpp
    title: cp-algo/algebra/common.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/algebra/modular.hpp
    title: cp-algo/algebra/modular.hpp
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
    PROBLEM: https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum
    document_title: Dynamic Range Affine Range Sum
    links:
    - https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum
  bundledCode: "#line 1 \"verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp\"\
    \n// @brief Dynamic Range Affine Range Sum\n#define PROBLEM \"https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum\"\
    \n#line 1 \"cp-algo/algebra/modular.hpp\"\n\n\n#line 1 \"cp-algo/random/rng.hpp\"\
    \n\n\n#include <chrono>\n#include <random>\nnamespace cp_algo::random {\n    std::mt19937_64\
    \ rng(std::chrono::steady_clock::now().time_since_epoch().count()); \n}\n\n#line\
    \ 1 \"cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n#include <cassert>\n\
    namespace cp_algo::algebra {\n    // a * x + b\n    template<typename base>\n\
    \    struct lin {\n        base a = 1, b = 0;\n        std::optional<base> c;\n\
    \        lin() {}\n        lin(base b): a(0), b(b) {}\n        lin(base a, base\
    \ b): a(a), b(b) {}\n        lin(base a, base b, base _c): a(a), b(b), c(_c) {}\n\
    \n        // polynomial product modulo x^2 - c\n        lin operator * (const\
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
    \ * A + b * B, c * A + d * B};\n        }\n    };\n}\n\n#line 1 \"cp-algo/algebra/common.hpp\"\
    \n\n\n#include <cstdint>\nnamespace cp_algo::algebra {\n    const int maxn = 1\
    \ << 20;\n    const int magic = 250; // threshold for sizes to run the naive algo\n\
    \n    auto bpow(auto x, int64_t n, auto ans) {\n        for(; n; n /= 2, x = x\
    \ * x) {\n            if(n % 2) {\n                ans = ans * x;\n          \
    \  }\n        }\n        return ans;\n    }\n    template<typename T>\n    T bpow(T\
    \ const& x, int64_t n) {\n        return bpow(x, n, T(1));\n    }\n\n    template<typename\
    \ T>\n    T fact(int n) {\n        static T F[maxn];\n        static bool init\
    \ = false;\n        if(!init) {\n            F[0] = T(1);\n            for(int\
    \ i = 1; i < maxn; i++) {\n                F[i] = F[i - 1] * T(i);\n         \
    \   }\n            init = true;\n        }\n        return F[n];\n    }\n    \n\
    \    template<typename T>\n    T rfact(int n) {\n        static T F[maxn];\n \
    \       static bool init = false;\n        if(!init) {\n            F[maxn - 1]\
    \ = T(1) / fact<T>(maxn - 1);\n            for(int i = maxn - 2; i >= 0; i--)\
    \ {\n                F[i] = F[i + 1] * T(i + 1);\n            }\n            init\
    \ = true;\n        }\n        return F[n];\n    }\n\n    template<typename T>\n\
    \    T small_inv(int n) {\n        static T F[maxn];\n        static bool init\
    \ = false;\n        if(!init) {\n            for(int i = 1; i < maxn; i++) {\n\
    \                F[i] = rfact<T>(i) * fact<T>(i - 1);\n            }\n       \
    \     init = true;\n        }\n        return F[n];\n    }\n}\n\n#line 6 \"cp-algo/algebra/modular.hpp\"\
    \n#include <algorithm>\n#include <iostream>\n#line 9 \"cp-algo/algebra/modular.hpp\"\
    \nnamespace cp_algo::algebra {\n    template<int m>\n    struct modular {\n  \
    \      // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm\n        std::optional<modular>\
    \ sqrt() const {\n            if(r == 0) {\n                return 0;\n      \
    \      } else if(bpow(*this, (m - 1) / 2) != modular(1)) {\n                return\
    \ std::nullopt;\n            } else {\n                while(true) {\n       \
    \             modular z = random::rng();\n                    if(z * z == *this)\
    \ {\n                        return z;\n                    }\n              \
    \      lin<modular> x(1, z, *this); // x + z (mod x^2 - b)\n                 \
    \   x = bpow(x, (m - 1) / 2, lin<modular>(0, 1, *this));\n                   \
    \ if(x.a != modular(0)) {\n                        return x.a.inv();\n       \
    \             }\n                }\n            }\n        }\n        \n     \
    \   uint64_t r;\n        constexpr modular(): r(0) {}\n        constexpr modular(int64_t\
    \ rr): r(rr % m) {r = std::min<uint64_t>(r, r + m);}\n        modular inv() const\
    \ {return bpow(*this, m - 2);}\n        modular operator - () const {return std::min(-r,\
    \ m - r);}\n        modular operator * (const modular &t) const {return r * t.r;}\n\
    \        modular operator / (const modular &t) const {return *this * t.inv();}\n\
    \        modular& operator += (const modular &t) {r += t.r; r = std::min<uint64_t>(r,\
    \ r - m); return *this;}\n        modular& operator -= (const modular &t) {r -=\
    \ t.r; r = std::min<uint64_t>(r, r + m); return *this;}\n        modular operator\
    \ + (const modular &t) const {return modular(*this) += t;}\n        modular operator\
    \ - (const modular &t) const {return modular(*this) -= t;}\n        modular& operator\
    \ *= (const modular &t) {return *this = *this * t;}\n        modular& operator\
    \ /= (const modular &t) {return *this = *this / t;}\n        \n        auto operator\
    \ <=> (const modular &t) const = default;\n        \n        explicit operator\
    \ int() const {return r;}\n        int64_t rem() const {return 2 * r > m ? r -\
    \ m : r;}\n\n        static constexpr uint64_t mm = (uint64_t)m * m;\n       \
    \ void add_unsafe(uint64_t t) {r += t; r = std::min<uint64_t>(r, r - mm);}\n \
    \       modular& normalize() {if(r >= m) r %= m; return *this;}\n    };\n    \n\
    \    template<int m>\n    std::istream& operator >> (std::istream &in, modular<m>\
    \ &x) {\n        return in >> x.r;\n    }\n    \n    template<int m>\n    std::ostream&\
    \ operator << (std::ostream &out, modular<m> const& x) {\n        return out <<\
    \ x.r % m;\n    }\n}\n\n#line 1 \"cp-algo/data_structures/treap/metas/reverse.hpp\"\
    \n\n\n#line 1 \"cp-algo/data_structures/treap/metas/base.hpp\"\n\n\n#line 1 \"\
    cp-algo/data_structures/treap/common.hpp\"\n\n\n#define _safe(t, op) (t ? t->op\
    \ : typename std::remove_reference_t<decltype(t->op)>())\n\n#line 4 \"cp-algo/data_structures/treap/metas/base.hpp\"\
    \n#include <functional>\n#line 7 \"cp-algo/data_structures/treap/metas/base.hpp\"\
    \n#define _safe_meta(i, op) _safe(i, _meta.op)\nnamespace cp_algo::data_structures::treap::metas\
    \ {\n    struct base_meta {\n        void pull(auto const, auto const){}\n   \
    \     void push(auto&, auto&){}\n    };\n}\n\n#line 6 \"cp-algo/data_structures/treap/metas/reverse.hpp\"\
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
    \n\n\n#line 5 \"cp-algo/data_structures/treap.hpp\"\n#include <array>\n/* Submissions\
    \ on Library Judge:\n  Range Reverse Range Sum, 558ms - https://judge.yosupo.jp/submission/147860\n\
    \  Cartesian Tree, 229ms - https://judge.yosupo.jp/submission/147858\n  Dynamic\
    \ Sequence Range Affine Range Sum, 2245ms - https://judge.yosupo.jp/submission/148948\n\
    */\nnamespace cp_algo::data_structures {\n    template<typename meta>\n    struct\
    \ treap_node {\n        using node = treap_node;\n        using treap = node*;\n\
    \        meta _meta;\n        int prior = random::rng();\n        size_t size\
    \ = 1;\n        treap children[2] = {nullptr, nullptr};\n        enum subtree\
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
    }\n\n#line 6 \"verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp\"\
    \n#include <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::data_structures;\n\
    \nusing base = cp_algo::algebra::modular<998244353>;\nusing meta = treap::metas::reverse_meta<base>;\n\
    using node = treap_node<meta>;\nusing treap_t = node::treap;\n\nvoid solve() {\n\
    \    istream_iterator<int> input(cin);\n    int n = *input++;\n    int q = *input++;\n\
    \    vector<treap_t> nodes(n);\n    generate_n(begin(nodes), n, [&](){\n     \
    \   return new node(meta(*input++));\n    });\n    auto me = node::build(nodes);\n\
    \n    while(q--) {\n        int t = *input++;\n        if(t == 0) {\n        \
    \    int i = *input++;\n            base x = *input++;\n            node::insert(me,\
    \ i, new node(meta(x)));\n        } else if(t == 1) {\n            node::erase(me,\
    \ *input++);\n        } else if(t == 2) {\n            int l = *input++;\n   \
    \         int r = *input++;\n            node::exec_on_segment(me, l, r, [](auto\
    \ &t) {\n                _safe_meta(t, reverse = 1);\n            });\n      \
    \  } else if(t == 3) {\n            int l = *input++;\n            int r = *input++;\n\
    \            base b = *input++;\n            base c = *input++;\n            node::exec_on_segment(me,\
    \ l, r, [b, c](auto &t) {\n                _safe_meta(t, add_push(meta::lin(b,\
    \ c)));\n            });\n        } else {\n            int l = *input++;\n  \
    \          int r = *input++;\n            node::exec_on_segment(me, l, r, [](auto\
    \ t) {\n                cout << _safe_meta(t, sum) << \"\\n\";\n            });\n\
    \        }\n    }\n}\n\nsigned main() {\n    //freopen(\"input.txt\", \"r\", stdin);\n\
    \    ios::sync_with_stdio(0);\n    cin.tie(0);\n    int t = 1;\n    while(t--)\
    \ {\n        solve();\n    }\n}\n"
  code: "// @brief Dynamic Range Affine Range Sum\n#define PROBLEM \"https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum\"\
    \n#include \"cp-algo/algebra/modular.hpp\"\n#include \"cp-algo/data_structures/treap/metas/reverse.hpp\"\
    \n#include \"cp-algo/data_structures/treap.hpp\"\n#include <bits/stdc++.h>\n\n\
    using namespace std;\nusing namespace cp_algo::data_structures;\n\nusing base\
    \ = cp_algo::algebra::modular<998244353>;\nusing meta = treap::metas::reverse_meta<base>;\n\
    using node = treap_node<meta>;\nusing treap_t = node::treap;\n\nvoid solve() {\n\
    \    istream_iterator<int> input(cin);\n    int n = *input++;\n    int q = *input++;\n\
    \    vector<treap_t> nodes(n);\n    generate_n(begin(nodes), n, [&](){\n     \
    \   return new node(meta(*input++));\n    });\n    auto me = node::build(nodes);\n\
    \n    while(q--) {\n        int t = *input++;\n        if(t == 0) {\n        \
    \    int i = *input++;\n            base x = *input++;\n            node::insert(me,\
    \ i, new node(meta(x)));\n        } else if(t == 1) {\n            node::erase(me,\
    \ *input++);\n        } else if(t == 2) {\n            int l = *input++;\n   \
    \         int r = *input++;\n            node::exec_on_segment(me, l, r, [](auto\
    \ &t) {\n                _safe_meta(t, reverse = 1);\n            });\n      \
    \  } else if(t == 3) {\n            int l = *input++;\n            int r = *input++;\n\
    \            base b = *input++;\n            base c = *input++;\n            node::exec_on_segment(me,\
    \ l, r, [b, c](auto &t) {\n                _safe_meta(t, add_push(meta::lin(b,\
    \ c)));\n            });\n        } else {\n            int l = *input++;\n  \
    \          int r = *input++;\n            node::exec_on_segment(me, l, r, [](auto\
    \ t) {\n                cout << _safe_meta(t, sum) << \"\\n\";\n            });\n\
    \        }\n    }\n}\n\nsigned main() {\n    //freopen(\"input.txt\", \"r\", stdin);\n\
    \    ios::sync_with_stdio(0);\n    cin.tie(0);\n    int t = 1;\n    while(t--)\
    \ {\n        solve();\n    }\n}\n"
  dependsOn:
  - cp-algo/algebra/modular.hpp
  - cp-algo/random/rng.hpp
  - cp-algo/algebra/affine.hpp
  - cp-algo/algebra/common.hpp
  - cp-algo/data_structures/treap/metas/reverse.hpp
  - cp-algo/data_structures/treap/metas/base.hpp
  - cp-algo/data_structures/treap/common.hpp
  - cp-algo/algebra/affine.hpp
  - cp-algo/data_structures/treap.hpp
  - cp-algo/random/rng.hpp
  - cp-algo/data_structures/treap/common.hpp
  isVerificationFile: true
  path: verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
  requiredBy: []
  timestamp: '2024-02-11 12:35:24+01:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
layout: document
redirect_from:
- /verify/verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp
- /verify/verify/data_structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp.html
title: Dynamic Range Affine Range Sum
---
