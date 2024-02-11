---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segtree.hpp
    title: cp-algo/data_structures/segtree.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segtree/metas/base.hpp
    title: cp-algo/data_structures/segtree/metas/base.hpp
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp
    title: cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum
    document_title: Range Chmin Chmax Add Range Sum
    links:
    - https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum
  bundledCode: "#line 1 \"verify/data_structures/segtree/range_chmin_chmax_add_range_sum.test.cpp\"\
    \n// @brief Range Chmin Chmax Add Range Sum\n#define PROBLEM \"https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum\"\
    \n#line 1 \"cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp\"\n\n\n\
    #line 1 \"cp-algo/data_structures/segtree/metas/base.hpp\"\n\n\nnamespace cp_algo::data_structures::segtree::metas\
    \ {\n    template<typename derived_meta>\n    struct base_meta {\n        using\
    \ meta = derived_meta;\n        virtual void pull(meta const&, meta const&, int,\
    \ int) {};\n        virtual void push(meta*, meta*, int, int) {};\n    };\n}\n\
    \n#line 4 \"cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp\"\n#include\
    \ <functional>\n#include <algorithm>\n#include <cstdint>\nnamespace cp_algo::data_structures::segtree::metas\
    \ {\n    struct chmin_chmax_sum_meta: base_meta<chmin_chmax_sum_meta> {\n    \
    \    static constexpr int64_t inf = 1e12;\n\n        using meta = chmin_chmax_sum_meta;\n\
    \        int64_t sum = 0, add = 0;\n\n        template<typename Comp>\n      \
    \  struct data {\n            int64_t val;\n            int64_t count = 1;\n \
    \           int64_t second = std::max(inf, -inf, comp);\n            static const\
    \ Comp comp;\n\n            data combine(data const& t) const {\n            \
    \    return comp(val, t.val) ? data{val, count, std::min(second, t.val, comp)}\n\
    \                        : comp(t.val, val) ? data{t.val, t.count, std::min(t.second,\
    \ val, comp)}\n                        : data{val, count + t.count, std::min(second,\
    \ t.second, comp)};\n            }\n\n            void add(int64_t b) {\n    \
    \            val += b;\n                second += b;\n            }\n\n      \
    \      int64_t normalize(int64_t L, int64_t R) {\n                int64_t old_val\
    \ = val;\n                val = std::clamp(val, L, R);\n                second\
    \ = std::clamp(second, L, R);\n                return count * (val - old_val);\n\
    \            }\n\n            bool stop(int64_t b) const {\n                return\
    \ !comp(val, b);\n            }\n            bool proceed(int64_t b) const {\n\
    \                return comp(b, second);\n            }\n        };\n        data<std::less<>>\
    \ mn = {sum};\n        data<std::greater<>> mx = {sum};\n        int64_t chmin\
    \ = inf, chmax = -inf;\n\n        chmin_chmax_sum_meta() {}\n        chmin_chmax_sum_meta(int64_t\
    \ val): sum(val) {}\n\n        void pull(meta const& L, meta const& R, int, int)\
    \ override {\n            sum = L.sum + R.sum;\n            mn = L.mn.combine(R.mn);\n\
    \            mx = L.mx.combine(R.mx);\n        }\n\n        void push(meta &t)\
    \ {\n            t.add += add; t.chmin += add; t.chmax += add;\n            t.chmin\
    \ = std::clamp(t.chmin, chmax, chmin);\n            t.chmax = std::clamp(t.chmax,\
    \ chmax, chmin);\n        }\n\n        void push(meta* L, meta* R, int l, int\
    \ r) override {\n            if(r - l > 1) {\n                push(*L);\n    \
    \            push(*R);\n            }\n            if(add) {\n               \
    \ sum += (r - l) * add;\n                mn.add(add);\n                mx.add(add);\n\
    \            }\n            bool same = mn.val == mx.val;\n            auto to_add\
    \ = mn.normalize(chmax, chmin) + mx.normalize(chmax, chmin);\n            sum\
    \ += same ? to_add / 2 : to_add;\n            if(mn.val == mx.val) {\n       \
    \         mx = {mx.val, r - l};\n                mn = {mn.val, r - l};\n     \
    \       }\n            add = 0;\n            chmin = inf;\n            chmax =\
    \ -inf;\n        }\n\n        static auto proceed_chmin(int64_t b) {\n       \
    \     return [b](meta const& t) {return t.mx.proceed(b);};\n        }\n      \
    \  static auto stop_chmin(int64_t b) {\n            return [b](meta const& t)\
    \ {return t.mx.stop(b);};\n        }\n        static auto proceed_chmax(int64_t\
    \ b) {\n            return [b](meta const& t) {return t.mn.proceed(b);};\n   \
    \     }\n        static auto stop_chmax(int64_t b) {\n            return [b](meta\
    \ const& t) {return t.mn.stop(b);};\n        }\n    };\n}\n\n#line 1 \"cp-algo/data_structures/segtree.hpp\"\
    \n\n\n#include <vector>\nnamespace cp_algo::data_structures {\n    template<typename\
    \ meta>\n    struct segtree_t {\n        const int N;\n        std::vector<meta>\
    \ _meta;\n\n        segtree_t(int n): N(n), _meta(4 * N) {}\n\n        segtree_t(std::vector<meta>\
    \ leafs): N(size(leafs)), _meta(4 * N) {\n            build(leafs);\n        }\n\
    \n        void pull(int v, int l, int r) {\n            if(r - l > 1) {\n    \
    \            _meta[v].pull(_meta[2 * v], _meta[2 * v + 1], l, r);\n          \
    \  }\n        }\n\n        void push(int v, int l, int r) {\n            if(r\
    \ - l > 1) {\n                _meta[v].push(&_meta[2 * v], &_meta[2 * v + 1],\
    \ l, r);\n            } else {\n                _meta[v].push(nullptr, nullptr,\
    \ l, r);\n            }\n        }\n\n        void build(auto &a, int v, size_t\
    \ l, size_t r) {\n            if(r - l == 1) {\n                if(l < size(a))\
    \ {\n                    _meta[v] = a[l];\n                }\n            } else\
    \ {\n                size_t m = (l + r) / 2;\n                build(a, 2 * v,\
    \ l, m);\n                build(a, 2 * v + 1, m, r);\n                pull(v,\
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
    \        }\n    };\n}\n\n#line 5 \"verify/data_structures/segtree/range_chmin_chmax_add_range_sum.test.cpp\"\
    \n#include <bits/stdc++.h>\n\nusing namespace std;\nusing namespace cp_algo::data_structures;\n\
    using meta = segtree::metas::chmin_chmax_sum_meta;\n\nvoid solve() {\n    int\
    \ n, q;\n    cin >> n >> q;\n    vector<meta> a(n);\n    for(int i = 0; i < n;\
    \ i++) {\n        int64_t ai;\n        cin >> ai;\n        a[i] = {ai};\n    }\n\
    \    segtree_t<meta> me(a);\n    while(q--) {\n        int t, l, r;\n        int64_t\
    \ b;\n        cin >> t >> l >> r;\n        if(t == 0) {\n            cin >> b;\n\
    \            me.exec_on_segment(l, r,\n                [b](auto& meta) {meta.chmin\
    \ = b;},\n                meta::proceed_chmin(b), meta::stop_chmin(b));\n    \
    \    } else if(t == 1) {\n            cin >> b;\n            me.exec_on_segment(l,\
    \ r,\n                [b](auto& meta) {meta.chmax = b;},\n                meta::proceed_chmax(b),\
    \ meta::stop_chmax(b));\n        } else if(t == 2) {\n            cin >> b;\n\
    \            me.exec_on_segment(l, r,\n                [b](auto& meta) {meta.add\
    \ = b;});\n        } else {\n            int64_t ans = 0;\n            me.exec_on_segment(l,\
    \ r, [&](auto& meta) {\n                ans += meta.sum;});\n            cout\
    \ << ans << \"\\n\";\n        }\n    }\n}\n\nsigned main() {\n    //freopen(\"\
    input.txt\", \"r\", stdin);\n    ios::sync_with_stdio(0);\n    cin.tie(0);\n \
    \   int t = 1;\n    while(t--) {\n        solve();\n    }\n}\n"
  code: "// @brief Range Chmin Chmax Add Range Sum\n#define PROBLEM \"https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum\"\
    \n#include \"cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp\"\n#include\
    \ \"cp-algo/data_structures/segtree.hpp\"\n#include <bits/stdc++.h>\n\nusing namespace\
    \ std;\nusing namespace cp_algo::data_structures;\nusing meta = segtree::metas::chmin_chmax_sum_meta;\n\
    \nvoid solve() {\n    int n, q;\n    cin >> n >> q;\n    vector<meta> a(n);\n\
    \    for(int i = 0; i < n; i++) {\n        int64_t ai;\n        cin >> ai;\n \
    \       a[i] = {ai};\n    }\n    segtree_t<meta> me(a);\n    while(q--) {\n  \
    \      int t, l, r;\n        int64_t b;\n        cin >> t >> l >> r;\n       \
    \ if(t == 0) {\n            cin >> b;\n            me.exec_on_segment(l, r,\n\
    \                [b](auto& meta) {meta.chmin = b;},\n                meta::proceed_chmin(b),\
    \ meta::stop_chmin(b));\n        } else if(t == 1) {\n            cin >> b;\n\
    \            me.exec_on_segment(l, r,\n                [b](auto& meta) {meta.chmax\
    \ = b;},\n                meta::proceed_chmax(b), meta::stop_chmax(b));\n    \
    \    } else if(t == 2) {\n            cin >> b;\n            me.exec_on_segment(l,\
    \ r,\n                [b](auto& meta) {meta.add = b;});\n        } else {\n  \
    \          int64_t ans = 0;\n            me.exec_on_segment(l, r, [&](auto& meta)\
    \ {\n                ans += meta.sum;});\n            cout << ans << \"\\n\";\n\
    \        }\n    }\n}\n\nsigned main() {\n    //freopen(\"input.txt\", \"r\", stdin);\n\
    \    ios::sync_with_stdio(0);\n    cin.tie(0);\n    int t = 1;\n    while(t--)\
    \ {\n        solve();\n    }\n}\n"
  dependsOn:
  - cp-algo/data_structures/segtree/metas/chmin_chmax_add.hpp
  - cp-algo/data_structures/segtree/metas/base.hpp
  - cp-algo/data_structures/segtree.hpp
  isVerificationFile: true
  path: verify/data_structures/segtree/range_chmin_chmax_add_range_sum.test.cpp
  requiredBy: []
  timestamp: '2024-02-11 15:34:32+01:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: verify/data_structures/segtree/range_chmin_chmax_add_range_sum.test.cpp
layout: document
redirect_from:
- /verify/verify/data_structures/segtree/range_chmin_chmax_add_range_sum.test.cpp
- /verify/verify/data_structures/segtree/range_chmin_chmax_add_range_sum.test.cpp.html
title: Range Chmin Chmax Add Range Sum
---
