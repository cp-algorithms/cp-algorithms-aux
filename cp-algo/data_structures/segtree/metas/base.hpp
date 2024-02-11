#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP
#define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP
namespace cp_algo::data_structures::segtree::metas {
    template<typename derived_meta>
    struct base_meta {
        using meta = derived_meta;
        virtual void pull(meta const&, meta const&, int, int) {};
        virtual void push(meta*, meta*, int, int) {};
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP