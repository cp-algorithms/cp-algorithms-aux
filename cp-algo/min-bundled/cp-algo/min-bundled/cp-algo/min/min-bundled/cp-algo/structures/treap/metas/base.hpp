#line 1 "cp-algo/min-bundled/cp-algo/min/min-bundled/cp-algo/structures/treap/metas/base.hpp"
#line 1 "cp-algo/min/min-bundled/cp-algo/structures/treap/metas/base.hpp"
#line 1 "cp-algo/structures/treap/metas/base.hpp"
#line 1 "cp-algo/structures/treap/common.hpp"
#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())
#line 4 "cp-algo/structures/treap/metas/base.hpp"
#include <functional>
#include <algorithm>
#include <cstdint>
#define _safe_meta(i, op) _safe(i, _meta.op)
namespace cp_algo::structures::treap::metas{struct base_meta{void pull(auto const,auto const){}
void push(auto&,auto&){}};}