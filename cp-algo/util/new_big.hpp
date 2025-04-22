#ifndef CP_ALGO_UTIL_NEW_BIG_HPP
#define CP_ALGO_UTIL_NEW_BIG_HPP
#include <sys/mman.h>
namespace cp_algo {
    template<typename T>
    auto new_big(size_t len) {
        auto raw = mmap(nullptr, len * sizeof(T),
            PROT_READ | PROT_WRITE,
            MAP_PRIVATE | MAP_ANONYMOUS,
            -1, 0);
        madvise(raw, len * sizeof(T), MADV_HUGEPAGE);
        return static_cast<T*>(raw);
    }
    template<typename T>
    void delete_big(T* ptr, size_t len) {
        munmap(ptr, len * sizeof(T));
    }
}
#endif // CP_ALGO_UTIL_NEW_BIG_HPP
