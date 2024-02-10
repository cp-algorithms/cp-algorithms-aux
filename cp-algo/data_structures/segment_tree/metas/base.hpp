#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP
#define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP
#include <functional>
#include <algorithm>
#include <cstdint>
namespace cp_algo::data_structures::segment_tree::metas {
    template<typename derived_meta>
    struct base_meta {
        using meta = derived_meta;
        void pull(meta const&, meta const&, int, int) {};
        void push(meta*, meta*, int, int) {};
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP