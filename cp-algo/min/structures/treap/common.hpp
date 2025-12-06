#ifndef CP_ALGO_STRUCTURES_TREAP_COMMON_HPP
#define CP_ALGO_STRUCTURES_TREAP_COMMON_HPP
#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())
#endif
