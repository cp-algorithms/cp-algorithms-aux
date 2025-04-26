#ifndef CP_ALGO_UTIL_big_alloc_HPP
#define CP_ALGO_UTIL_big_alloc_HPP

#include <cstddef>
#include <iostream>

// Single macro to detect POSIX platforms (Linux, Unix, macOS)
#if defined(__linux__) || defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
#endif

namespace cp_algo {
    template <typename T>
    class big_alloc: public std::allocator<T> {
        public:
        using value_type = T;
        using base = std::allocator<T>;

        big_alloc() noexcept = default;

        template <typename U>
        big_alloc(const big_alloc<U>&) noexcept {}

#if CP_ALGO_USE_MMAP
        [[nodiscard]] T* allocate(std::size_t n) {
            if(n * sizeof(T) < 4 * 1024 * 1024) {
                return base::allocate(n);
            }
            n *= sizeof(T);
            void* raw = mmap(nullptr, n,
                            PROT_READ | PROT_WRITE,
                            MAP_PRIVATE | MAP_ANONYMOUS,
                            -1, 0);
            madvise(raw, n, MADV_HUGEPAGE);
            madvise(raw, n, MADV_POPULATE_WRITE);
            return static_cast<T*>(raw);
        }
#endif

#if CP_ALGO_USE_MMAP
        void deallocate(T* p, std::size_t n) noexcept {
            if(n * sizeof(T) < 4 * 1024 * 1024) {
                return base::deallocate(p, n);
            }
            if(p) {
                munmap(p, n * sizeof(T));
            }
        }
#endif
    };
}
#endif // CP_ALGO_UTIL_big_alloc_HPP
