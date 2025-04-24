#ifndef CP_ALGO_UTIL_NEW_BIG_HPP
#define CP_ALGO_UTIL_NEW_BIG_HPP
#include <cstddef>

// Single macro to detect POSIX platforms (Linux, Unix, macOS)
#if defined(__linux__) || defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
#endif

namespace cp_algo {
    template<typename T>
    T* new_big(std::size_t len) {
#if CP_ALGO_USE_MMAP
        void* raw = mmap(nullptr, len * sizeof(T),
                         PROT_READ | PROT_WRITE,
                         MAP_PRIVATE | MAP_ANONYMOUS,
                         -1, 0);
        // Advise the kernel for huge pages and pre-populate writes
        madvise(raw, len * sizeof(T), MADV_HUGEPAGE);
        madvise(raw, len * sizeof(T), MADV_POPULATE_WRITE);
        return static_cast<T*>(raw);
#else
        // Fallback allocation for non-POSIX platforms
        return new T[len];
#endif
    }
    template<typename T>
    void delete_big(T* ptr, [[maybe_unused]] std::size_t len) {
#if CP_ALGO_USE_MMAP
        munmap(ptr, len * sizeof(T));
#else
        // Match allocation with delete[] on non-POSIX platforms
        delete[] ptr;
#endif
    }
}
#endif // CP_ALGO_UTIL_NEW_BIG_HPP
