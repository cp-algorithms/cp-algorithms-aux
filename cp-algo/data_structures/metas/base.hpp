#ifndef DATA_STRUCTURES_METAS_BASE_HPP
#define DATA_STRUCTURES_METAS_BASE_HPP
#include <functional>
#include <algorithm>
#include <cstdint>
namespace data_structures {
    namespace segment_tree {
        namespace metas {
            template<typename derived_meta>
            struct base_meta {
                using meta = derived_meta;
                virtual void pull(meta const&, meta const&, int, int) = 0;
                virtual void push(meta*, meta*, int, int) = 0;
            };
        }
    }
}
#endif // DATA_STRUCTURES_METAS_BASE_HPP