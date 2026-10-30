// Licensed under GPL-3.0-only with the online-judge exception.
// See cp-algo/LICENSE and cp-algo/ONLINE-JUDGE-EXCEPTION.
#ifndef CP_ALGO_UTIL_BUMP_ALLOC_HPP
#define CP_ALGO_UTIL_BUMP_ALLOC_HPP
#include <cstddef>
#include "big_alloc.hpp"
namespace cp_algo {
    template<class T, size_t max_len>
    struct bump_alloc {
        static inline char* buf = big_alloc<char>().allocate(max_len * sizeof(T));
        static inline size_t buf_ind = max_len * sizeof(T);
        using value_type = T;
        template <class U> struct rebind { using other = bump_alloc<U, max_len>; };
        constexpr bool operator==(const bump_alloc&) const = default;
        bump_alloc() = default;
        template<class U> bump_alloc(const U&) {}
        T* allocate(size_t n) {
            buf_ind -= n * sizeof(T);
            buf_ind &= 0 - alignof(T);
            return (T*)(buf + buf_ind);
        }
        void deallocate(T*, size_t) {}
    };
}
#endif // CP_ALGO_UTIL_BUMP_ALLOC_HPP
