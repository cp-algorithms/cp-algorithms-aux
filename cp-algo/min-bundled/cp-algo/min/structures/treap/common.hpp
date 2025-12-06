#line 1 "cp-algo/min/structures/treap/common.hpp"
#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())