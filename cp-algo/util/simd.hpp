#ifndef CP_ALGO_UTIL_SIMD_HPP
#define CP_ALGO_UTIL_SIMD_HPP
#include <experimental/simd>
#include <cstdint>
#include <cstddef>
#include <memory>
namespace cp_algo {
    template<typename T, size_t len>
    using simd [[gnu::vector_size(len * sizeof(T))]] = T;
    using u32x8 = simd<uint32_t, 8>;
    using i64x4 = simd<int64_t, 4>;
    using u64x4 = simd<uint64_t, 4>;
    using i32x4 = simd<int32_t, 4>;
    using u32x4 = simd<uint32_t, 4>;
    using i16x4 = simd<int16_t, 4>;
    using dx4 = simd<double, 4>;

    [[gnu::always_inline]] inline dx4 abs(dx4 a) {
    return a < 0 ? -a : a;
    }

    // https://stackoverflow.com/a/77376595
    // works for ints in (-2^51, 2^51)
    static constexpr dx4 magic = dx4() + (3ULL << 51);
    [[gnu::always_inline]] inline i64x4 lround(dx4 x) {
        return i64x4(x + magic) - i64x4(magic);
    }
    [[gnu::always_inline]] inline dx4 to_double(i64x4 x) {
        return dx4(x + i64x4(magic)) - magic;
    }

    [[gnu::always_inline]] inline dx4 round(dx4 a) {
        return dx4{
            std::nearbyint(a[0]),
            std::nearbyint(a[1]),
            std::nearbyint(a[2]),
            std::nearbyint(a[3])
        };
    }

    [[gnu::always_inline]] inline u64x4 low32(u64x4 x) {
        return x & uint32_t(-1);
    }

    [[gnu::always_inline]] inline u64x4 montgomery_reduce(u64x4 x, uint32_t mod, uint32_t imod) {
        auto x_ninv = u64x4(u32x8(x) * (u32x8() + imod));
#ifdef __AVX2__
        x += u64x4(_mm256_mul_epu32(__m256i(x_ninv), __m256i() + mod));
#else
        x += low32(x_ninv) * mod;
#endif
        return x >> 32;
    }

    [[gnu::always_inline]] inline u64x4 montgomery_mul(u64x4 x, u64x4 y, uint32_t mod, uint32_t imod) {
#ifdef __AVX2__
        return montgomery_reduce(u64x4(_mm256_mul_epu32(__m256i(x), __m256i(y))), mod, imod);
#else
        return montgomery_reduce(low32(x) * low32(y), mod, imod);
#endif
    }

    [[gnu::always_inline]] inline u32x8 montgomery_mul(u32x8 x, u32x8 y, uint32_t mod, uint32_t imod) {
        auto x0246 = u64x4(x);
        auto y0246 = u64x4(y);
        auto x1357 = u64x4(x) >> 32;
        auto y1357 = u64x4(y) >> 32;
        return u32x8(montgomery_mul(x0246, y0246, mod, imod)) |
               u32x8(montgomery_mul(x1357, y1357, mod, imod) << 32);
    }

    [[gnu::always_inline]] inline dx4 rotate_right(dx4 x) {
        static constexpr u64x4 shuffler = {3, 0, 1, 2};
        return __builtin_shuffle(x, shuffler);
    }

    template<std::size_t Align = 32>
    [[gnu::always_inline]] inline bool is_aligned(const auto* p) noexcept {
        return (reinterpret_cast<std::uintptr_t>(p) % Align) == 0;
    }

    template<class Target>
    [[gnu::always_inline]] inline Target& vector_cast(auto &&p) {
        return *reinterpret_cast<Target*>(std::assume_aligned<alignof(Target)>(&p));
    }
}
#endif // CP_ALGO_UTIL_SIMD_HPP
