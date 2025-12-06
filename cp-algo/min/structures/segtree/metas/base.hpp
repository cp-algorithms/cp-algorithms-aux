#ifndef CP_ALGO_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP
#define CP_ALGO_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP
#include <cstddef>
namespace cp_algo::structures::segtree::metas{
template<typename derived_meta>
struct base_meta{
using meta=derived_meta;
virtual void pull(meta const&,meta const&,size_t,size_t){};
virtual void push(meta*,meta*,size_t,size_t){};
};
}
#endif
