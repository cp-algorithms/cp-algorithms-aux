#line 1 "cp-algo/min-bundled/cp-algo/structures/treap/common.hpp"
#line 1 "cp-algo/structures/treap/common.hpp"
#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())