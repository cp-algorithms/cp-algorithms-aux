---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: cp-algo/data_structures/segment_tree/metas/base.hpp
    title: cp-algo/data_structures/segment_tree/metas/base.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: verify/data_structures/segment_tree/range_chmin_chmax_add_range_sum.test.cpp
    title: Range Chmin Chmax Add Range Sum
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"cp-algo/data_structures/segment_tree/metas/chmin_chmax_add.hpp\"\
    \n\n\n#line 1 \"cp-algo/data_structures/segment_tree/metas/base.hpp\"\n\n\n#include\
    \ <functional>\n#include <algorithm>\n#include <cstdint>\nnamespace cp_algo::data_structures::segment_tree::metas\
    \ {\n    template<typename derived_meta>\n    struct base_meta {\n        using\
    \ meta = derived_meta;\n        virtual void pull(meta const&, meta const&, int,\
    \ int) = 0;\n        virtual void push(meta*, meta*, int, int) = 0;\n    };\n\
    }\n\n#line 7 \"cp-algo/data_structures/segment_tree/metas/chmin_chmax_add.hpp\"\
    \nnamespace cp_algo::data_structures::segment_tree::metas {\n    struct chmin_chmax_sum_meta:\
    \ base_meta<chmin_chmax_sum_meta> {\n        static constexpr int64_t inf = 1e12;\n\
    \n        using meta = chmin_chmax_sum_meta;\n        int64_t sum = 0, add = 0;\n\
    \n        template<typename Comp>\n        struct data {\n            int64_t\
    \ val;\n            int64_t count = 1;\n            int64_t second = std::max(inf,\
    \ -inf, comp);\n            static const Comp comp;\n\n            data combine(data\
    \ const& t) const {\n                return comp(val, t.val) ? data{val, count,\
    \ std::min(second, t.val, comp)}\n                        : comp(t.val, val) ?\
    \ data{t.val, t.count, std::min(t.second, val, comp)}\n                      \
    \  : data{val, count + t.count, std::min(second, t.second, comp)};\n         \
    \   }\n\n            void add(int64_t b) {\n                val += b;\n      \
    \          second += b;\n            }\n\n            int64_t normalize(int64_t\
    \ L, int64_t R) {\n                int64_t old_val = val;\n                val\
    \ = std::clamp(val, L, R);\n                second = std::clamp(second, L, R);\n\
    \                return count * (val - old_val);\n            }\n\n          \
    \  bool stop(int64_t b) const {\n                return !comp(val, b);\n     \
    \       }\n            bool proceed(int64_t b) const {\n                return\
    \ comp(b, second);\n            }\n        };\n        data<std::less<int64_t>>\
    \ mn = {sum};\n        data<std::greater<int64_t>> mx = {sum};\n        int64_t\
    \ chmin = inf, chmax = -inf;\n\n        chmin_chmax_sum_meta() {}\n        chmin_chmax_sum_meta(int64_t\
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
    \ const& t) {return t.mn.stop(b);};\n        }\n    };\n}\n\n"
  code: "#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_CHMIN_CHMAX_ADD_HPP\n\
    #define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_CHMIN_CHMAX_ADD_HPP\n#include\
    \ \"base.hpp\"\n#include <functional>\n#include <algorithm>\n#include <cstdint>\n\
    namespace cp_algo::data_structures::segment_tree::metas {\n    struct chmin_chmax_sum_meta:\
    \ base_meta<chmin_chmax_sum_meta> {\n        static constexpr int64_t inf = 1e12;\n\
    \n        using meta = chmin_chmax_sum_meta;\n        int64_t sum = 0, add = 0;\n\
    \n        template<typename Comp>\n        struct data {\n            int64_t\
    \ val;\n            int64_t count = 1;\n            int64_t second = std::max(inf,\
    \ -inf, comp);\n            static const Comp comp;\n\n            data combine(data\
    \ const& t) const {\n                return comp(val, t.val) ? data{val, count,\
    \ std::min(second, t.val, comp)}\n                        : comp(t.val, val) ?\
    \ data{t.val, t.count, std::min(t.second, val, comp)}\n                      \
    \  : data{val, count + t.count, std::min(second, t.second, comp)};\n         \
    \   }\n\n            void add(int64_t b) {\n                val += b;\n      \
    \          second += b;\n            }\n\n            int64_t normalize(int64_t\
    \ L, int64_t R) {\n                int64_t old_val = val;\n                val\
    \ = std::clamp(val, L, R);\n                second = std::clamp(second, L, R);\n\
    \                return count * (val - old_val);\n            }\n\n          \
    \  bool stop(int64_t b) const {\n                return !comp(val, b);\n     \
    \       }\n            bool proceed(int64_t b) const {\n                return\
    \ comp(b, second);\n            }\n        };\n        data<std::less<int64_t>>\
    \ mn = {sum};\n        data<std::greater<int64_t>> mx = {sum};\n        int64_t\
    \ chmin = inf, chmax = -inf;\n\n        chmin_chmax_sum_meta() {}\n        chmin_chmax_sum_meta(int64_t\
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
    \ const& t) {return t.mn.stop(b);};\n        }\n    };\n}\n#endif // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_CHMIN_CHMAX_ADD_HPP"
  dependsOn:
  - cp-algo/data_structures/segment_tree/metas/base.hpp
  isVerificationFile: false
  path: cp-algo/data_structures/segment_tree/metas/chmin_chmax_add.hpp
  requiredBy: []
  timestamp: '2024-02-10 22:49:03+01:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - verify/data_structures/segment_tree/range_chmin_chmax_add_range_sum.test.cpp
documentation_of: cp-algo/data_structures/segment_tree/metas/chmin_chmax_add.hpp
layout: document
redirect_from:
- /library/cp-algo/data_structures/segment_tree/metas/chmin_chmax_add.hpp
- /library/cp-algo/data_structures/segment_tree/metas/chmin_chmax_add.hpp.html
title: cp-algo/data_structures/segment_tree/metas/chmin_chmax_add.hpp
---
