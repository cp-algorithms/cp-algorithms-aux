#ifndef CP_ALGO_DATA_STRUCTURES_TREAP_METAS_BASE_HPP
#define CP_ALGO_DATA_STRUCTURES_TREAP_METAS_BASE_HPP
#include "../../treap.hpp"
#include <functional>
#include <algorithm>
#include <cstdint>
#define _safe_meta(i, op) _safe(i, _meta.op)
namespace cp_algo::data_structures::treap::metas {
    struct base_meta {
        void pull(auto const, auto const);
        void push(auto&, auto&);
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_TREAP_METAS_BASE_HPP