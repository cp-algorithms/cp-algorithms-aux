#ifndef CP_ALGO_UTIL_BUMP_ALLOC_HPP
#define CP_ALGO_UTIL_BUMP_ALLOC_HPP
#include <cstddef>
namespace cp_algo {
    char buf[450 << 20] alignas(32);
    size_t buf_ind = sizeof buf;
    template<class T> struct bump_alloc {
        typedef T value_type;
        bump_alloc() {}
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
