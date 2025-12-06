#line 1 "cp-algo/min-bundled/cp-algo/min/min-bundled/cp-algo/structures/segtree/metas/base.hpp"
#line 1 "cp-algo/min/min-bundled/cp-algo/structures/segtree/metas/base.hpp"
#line 1 "cp-algo/structures/segtree/metas/base.hpp"
#include <cstddef>
namespace cp_algo::structures::segtree::metas{template<typename derived_meta>
struct base_meta{using meta=derived_meta;
virtual void pull(meta const&,meta const&,size_t,size_t){};
virtual void push(meta*,meta*,size_t,size_t){};};}