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
    path: cp-algo/data_structures/segment_tree.hpp
    title: cp-algo/data_structures/segment_tree.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segment_tree/metas/affine.hpp
    title: cp-algo/data_structures/segment_tree/metas/affine.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segment_tree/metas/base.hpp
    title: cp-algo/data_structures/segment_tree/metas/base.hpp
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
    PROBLEM: https://judge.yosupo.jp/problem/range_affine_range_sum
    document_title: Range Affine Range Sum
    links:
    - https://judge.yosupo.jp/problem/range_affine_range_sum
  bundledCode: "#line 1 \"verify/data_structures/segment_tree/range_affine_range_sum.test.cpp\"\
    \n// @brief Range Affine Range Sum\n#define PROBLEM \"https://judge.yosupo.jp/problem/range_affine_range_sum\"\
    \n#line 1 \"cp-algo/data_structures/segment_tree/metas/affine.hpp\"\n\n\n#line\
    \ 1 \"cp-algo/data_structures/segment_tree/metas/base.hpp\"\n\n\n#include <functional>\n\
    #include <algorithm>\n#include <cstdint>\nnamespace cp_algo::data_structures::segment_tree::metas\
    \ {\n    template<typename derived_meta>\n    struct base_meta {\n        using\
    \ meta = derived_meta;\n        virtual void pull(meta const&, meta const&, int,\
    \ int) {};\n        virtual void push(meta*, meta*, int, int) {};\n    };\n}\n\
    \n#line 1 \"cp-algo/algebra/affine.hpp\"\n\n\n#include <optional>\n#include <cassert>\n\
    namespace cp_algo::algebra {\n    template<typename base>\n    // a * x + b\n\
    \    struct lin {\n        base a = 1, b = 0;\n        std::optional<base> c;\n\
    \        lin() {}\n        lin(base b): a(0), b(b) {}\n        lin(base a, base\
    \ b): a(a), b(b) {}\n        lin(base a, base b, base _c): a(a), b(b), c(_c) {}\n\
    \n        // polynomial product modulo x^2 - c\n        lin operator * (const\
    \ lin& t) {\n            assert(c && t.c && *c == *t.c);\n            return {a\
    \ * t.b + b * t.a, b * t.b + a * t.a * (*c), *c};\n        }\n\n        // a *\
    \ (t.a * x + t.b) + b\n        lin apply(lin const& t) const {\n            return\
    \ {a * t.a, a * t.b + b};\n        }\n\n        void prepend(lin const& t) {\n\
    \            *this = t.apply(*this);\n        }\n\n        base eval(base x) const\
    \ {\n            return a * x + b;\n        }\n    };\n}\n\n#line 5 \"cp-algo/data_structures/segment_tree/metas/affine.hpp\"\
    \nnamespace cp_algo::data_structures::segment_tree::metas {\n    template<typename\
    \ base>\n    struct affine_meta: base_meta<affine_meta<base>> {\n        using\
    \ meta = affine_meta;\n        using lin = algebra::lin<base>;\n\n        base\
    \ sum = 0;\n        lin to_push = {};\n\n        affine_meta() {}\n        affine_meta(base\
    \ sum): sum(sum) {}\n\n        void push(meta *L, meta *R, int l, int r) override\
    \ {\n            if(to_push.a != 1 || to_push.b != 0) {\n                sum =\
    \ to_push.a * sum + to_push.b * (r - l);\n                if(r - l > 1) {\n  \
    \                  L->to_push.prepend(to_push);\n                    R->to_push.prepend(to_push);\n\
    \                }\n                to_push = {};\n            }\n        }\n\n\
    \        void pull(meta const& L, meta const& R, int, int) override {\n      \
    \      sum = L.sum + R.sum;\n        }\n    };\n}\n\n#line 1 \"cp-algo/data_structures/segment_tree.hpp\"\
    \n\n\n#include <vector>\nnamespace cp_algo::data_structures::segment_tree {\n\
    \    template<typename meta>\n    struct segment_tree {\n        const int N;\n\
    \        std::vector<meta> _meta;\n\n        segment_tree(int n): N(n), _meta(4\
    \ * N) {}\n\n        segment_tree(std::vector<meta> leafs): N(size(leafs)), _meta(4\
    \ * N) {\n            build(leafs);\n        }\n\n        void pull(int v, int\
    \ l, int r) {\n            if(r - l > 1) {\n                _meta[v].pull(_meta[2\
    \ * v], _meta[2 * v + 1], l, r);\n            }\n        }\n\n        void push(int\
    \ v, int l, int r) {\n            if(r - l > 1) {\n                _meta[v].push(&_meta[2\
    \ * v], &_meta[2 * v + 1], l, r);\n            } else {\n                _meta[v].push(nullptr,\
    \ nullptr, l, r);\n            }\n        }\n\n        void build(auto &a, int\
    \ v, size_t l, size_t r) {\n            if(r - l == 1) {\n                if(l\
    \ < size(a)) {\n                    _meta[v] = a[l];\n                }\n    \
    \        } else {\n                size_t m = (l + r) / 2;\n                build(a,\
    \ 2 * v, l, m);\n                build(a, 2 * v + 1, m, r);\n                pull(v,\
    \ l, r);\n            }\n        }\n\n        void build(auto &a) {\n        \
    \    build(a, 1, 0, N);\n        }\n\n        void exec_on_segment(int a, int\
    \ b, auto func, auto proceed, auto stop, int v, int l, int r) {\n            push(v,\
    \ l, r);\n            if(r <= a || b <= l || stop(_meta[v])) {\n             \
    \   return;\n            } else if(a <= l && r <= b && proceed(_meta[v])) {\n\
    \                func(_meta[v]);\n                push(v, l, r);\n           \
    \ } else {\n                int m = (l + r) / 2;\n                exec_on_segment(a,\
    \ b, func, proceed, stop, 2 * v, l, m);\n                exec_on_segment(a, b,\
    \ func, proceed, stop, 2 * v + 1, m, r);\n                pull(v, l, r);\n   \
    \         }\n        }\n\n        static constexpr auto default_true = [](auto\
    \ const&){return true;};\n        static constexpr auto default_false = [](auto\
    \ const&){return false;};\n\n        void exec_on_segment(int a, int b, auto func,\
    \ auto proceed, auto stop) {\n            exec_on_segment(a, b, func, proceed,\
    \ stop, 1, 0, N);\n        }\n\n        void exec_on_segment(int a, int b, auto\
    \ func) {\n            exec_on_segment(a, b, func, default_true, default_false);\n\
    \        }\n    };\n}\n\n#line 1 \"cp-algo/algebra/modular.hpp\"\n\n\n#line 1\
    \ \"cp-algo/random/rng.hpp\"\n\n\n#include <chrono>\n#include <random>\nnamespace\
    \ cp_algo::random {\n    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());\
    \ \n}\n\n#line 1 \"cp-algo/algebra/common.hpp\"\n\n\n#line 4 \"cp-algo/algebra/common.hpp\"\
    \nnamespace cp_algo::algebra {\n    const int maxn = 1 << 20;\n    const int magic\
    \ = 250; // threshold for sizes to run the naive algo\n\n    auto bpow(auto x,\
    \ int64_t n, auto ans) {\n        for(; n; n /= 2, x = x * x) {\n            if(n\
    \ % 2) {\n                ans = ans * x;\n            }\n        }\n        return\
    \ ans;\n    }\n    template<typename T>\n    T bpow(T const& x, int64_t n) {\n\
    \        return bpow(x, n, T(1));\n    }\n\n    template<typename T>\n    T fact(int\
    \ n) {\n        static T F[maxn];\n        static bool init = false;\n       \
    \ if(!init) {\n            F[0] = T(1);\n            for(int i = 1; i < maxn;\
    \ i++) {\n                F[i] = F[i - 1] * T(i);\n            }\n           \
    \ init = true;\n        }\n        return F[n];\n    }\n    \n    template<typename\
    \ T>\n    T rfact(int n) {\n        static T F[maxn];\n        static bool init\
    \ = false;\n        if(!init) {\n            F[maxn - 1] = T(1) / fact<T>(maxn\
    \ - 1);\n            for(int i = maxn - 2; i >= 0; i--) {\n                F[i]\
    \ = F[i + 1] * T(i + 1);\n            }\n            init = true;\n        }\n\
    \        return F[n];\n    }\n\n    template<typename T>\n    T small_inv(int\
    \ n) {\n        static T F[maxn];\n        static bool init = false;\n       \
    \ if(!init) {\n            for(int i = 1; i < maxn; i++) {\n                F[i]\
    \ = rfact<T>(i) * fact<T>(i - 1);\n            }\n            init = true;\n \
    \       }\n        return F[n];\n    }\n}\n\n#line 7 \"cp-algo/algebra/modular.hpp\"\
    \n#include <iostream>\n#line 9 \"cp-algo/algebra/modular.hpp\"\nnamespace cp_algo::algebra\
    \ {\n    template<int m>\n    struct modular {\n        // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm\n\
    \        std::optional<modular> sqrt() const {\n            if(r == 0) {\n   \
    \             return 0;\n            } else if(bpow(*this, (m - 1) / 2) != modular(1))\
    \ {\n                return std::nullopt;\n            } else {\n            \
    \    while(true) {\n                    modular z = random::rng();\n         \
    \           if(z * z == *this) {\n                        return z;\n        \
    \            }\n                    lin<modular> x(1, z, *this); // x + z (mod\
    \ x^2 - b)\n                    x = bpow(x, (m - 1) / 2, lin<modular>(0, 1, *this));\n\
    \                    if(x.a != modular(0)) {\n                        return x.a.inv();\n\
    \                    }\n                }\n            }\n        }\n        \n\
    \        uint64_t r;\n        constexpr modular(): r(0) {}\n        constexpr\
    \ modular(int64_t rr): r(rr % m) {r = std::min<uint64_t>(r, r + m);}\n       \
    \ modular inv() const {return bpow(*this, m - 2);}\n        modular operator -\
    \ () const {return std::min(-r, m - r);}\n        modular operator * (const modular\
    \ &t) const {return r * t.r;}\n        modular operator / (const modular &t) const\
    \ {return *this * t.inv();}\n        modular& operator += (const modular &t) {r\
    \ += t.r; r = std::min<uint64_t>(r, r - m); return *this;}\n        modular& operator\
    \ -= (const modular &t) {r -= t.r; r = std::min<uint64_t>(r, r + m); return *this;}\n\
    \        modular operator + (const modular &t) const {return modular(*this) +=\
    \ t;}\n        modular operator - (const modular &t) const {return modular(*this)\
    \ -= t;}\n        modular& operator *= (const modular &t) {return *this = *this\
    \ * t;}\n        modular& operator /= (const modular &t) {return *this = *this\
    \ / t;}\n        \n        auto operator <=> (const modular &t) const = default;\n\
    \        \n        explicit operator int() const {return r;}\n        int64_t\
    \ rem() const {return 2 * r > m ? r - m : r;}\n\n        static constexpr uint64_t\
    \ mm = (uint64_t)m * m;\n        void add_unsafe(uint64_t t) {r += t; r = std::min<uint64_t>(r,\
    \ r - mm);}\n        modular& normalize() {if(r >= m) r %= m; return *this;}\n\
    \    };\n    \n    template<int m>\n    std::istream& operator >> (std::istream\
    \ &in, modular<m> &x) {\n        return in >> x.r;\n    }\n    \n    template<int\
    \ m>\n    std::ostream& operator << (std::ostream &out, modular<m> const& x) {\n\
    \        return out << x.r % m;\n    }\n}\n\n#line 6 \"verify/data_structures/segment_tree/range_affine_range_sum.test.cpp\"\
    \n#include <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::data_structures::segment_tree;\n\
    \nusing base = cp_algo::algebra::modular<998244353>;\nusing meta = metas::affine_meta<base>;\n\
    \nvoid solve() {\n    int n, q;\n    cin >> n >> q;\n    vector<meta> a(n);\n\
    \    for(int i = 0; i < n; i++) {\n        int ai;\n        cin >> ai;\n     \
    \   a[i] = {ai};\n    }\n    segment_tree<meta> me(a);\n    while(q--) {\n   \
    \     int t;\n        cin >> t;\n        if(t == 0) {\n            int l, r, b,\
    \ c;\n            cin >> l >> r >> b >> c;\n            me.exec_on_segment(l,\
    \ r, [&](auto& meta) {\n                meta.to_push.prepend(meta::lin(b, c));\n\
    \            });\n        } else {\n            int l, r;\n            cin >>\
    \ l >> r;\n            base ans = 0;\n            me.exec_on_segment(l, r, [&](auto\
    \ meta) {\n                ans += meta.sum;\n            });\n            cout\
    \ << ans << \"\\n\";\n        }\n    }\n}\n\nsigned main() {\n    //freopen(\"\
    input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n    cin.tie(0);\n \
    \   int t = 1;\n    while(t--) {\n        solve();\n    }\n}\n"
  code: "// @brief Range Affine Range Sum\n#define PROBLEM \"https://judge.yosupo.jp/problem/range_affine_range_sum\"\
    \n#include \"cp-algo/data_structures/segment_tree/metas/affine.hpp\"\n#include\
    \ \"cp-algo/data_structures/segment_tree.hpp\"\n#include \"cp-algo/algebra/modular.hpp\"\
    \n#include <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::data_structures::segment_tree;\n\
    \nusing base = cp_algo::algebra::modular<998244353>;\nusing meta = metas::affine_meta<base>;\n\
    \nvoid solve() {\n    int n, q;\n    cin >> n >> q;\n    vector<meta> a(n);\n\
    \    for(int i = 0; i < n; i++) {\n        int ai;\n        cin >> ai;\n     \
    \   a[i] = {ai};\n    }\n    segment_tree<meta> me(a);\n    while(q--) {\n   \
    \     int t;\n        cin >> t;\n        if(t == 0) {\n            int l, r, b,\
    \ c;\n            cin >> l >> r >> b >> c;\n            me.exec_on_segment(l,\
    \ r, [&](auto& meta) {\n                meta.to_push.prepend(meta::lin(b, c));\n\
    \            });\n        } else {\n            int l, r;\n            cin >>\
    \ l >> r;\n            base ans = 0;\n            me.exec_on_segment(l, r, [&](auto\
    \ meta) {\n                ans += meta.sum;\n            });\n            cout\
    \ << ans << \"\\n\";\n        }\n    }\n}\n\nsigned main() {\n    //freopen(\"\
    input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n    cin.tie(0);\n \
    \   int t = 1;\n    while(t--) {\n        solve();\n    }\n}\n"
  dependsOn:
  - cp-algo/data_structures/segment_tree/metas/affine.hpp
  - cp-algo/data_structures/segment_tree/metas/base.hpp
  - cp-algo/algebra/affine.hpp
  - cp-algo/data_structures/segment_tree.hpp
  - cp-algo/algebra/modular.hpp
  - cp-algo/random/rng.hpp
  - cp-algo/algebra/affine.hpp
  - cp-algo/algebra/common.hpp
  isVerificationFile: true
  path: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
  requiredBy: []
  timestamp: '2024-02-11 00:23:03+01:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
layout: document
redirect_from:
- /verify/verify/data_structures/segment_tree/range_affine_range_sum.test.cpp
- /verify/verify/data_structures/segment_tree/range_affine_range_sum.test.cpp.html
title: Range Affine Range Sum
---
