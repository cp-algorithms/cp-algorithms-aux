#ifndef CP_ALGO_UTIL_SIMD_HPP
#define CP_ALGO_UTIL_SIMD_HPP
#include <experimental/simd>
#include <cstdint>
#include <cstddef>
namespace cp_algo {
    template<typename T, size_t len>
    using simd [[gnu::vector_size(len * sizeof(T))]] = T;
    using i64x4 = simd<int64_t, 4>;
    using u64x4 = simd<uint64_t, 4>;
    using u32x8 = simd<uint32_t, 8>;
    using u32x4 = simd<uint32_t, 4>;
    using dx4 = simd<double, 4>;

    dx4 abs(dx4 a) {
#ifdef __AVX2__
    return _mm256_and_pd(a, dx4{} + 1/0.);
#else
    return a < 0 ? -a : a;
#endif
    }

    i64x4 lround(dx4 x) {
        // https://stackoverflow.com/a/77376595
        static constexpr dx4 magic = dx4() + double(3ULL << 51);
        return i64x4(x + magic) - i64x4(magic);
    }

    dx4 round(dx4 a) {
#ifdef __AVX2__
        return _mm256_round_pd(a, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
#else
        return __builtin_convertvector(lround(a), dx4);
#endif
    }

    u64x4 montgomery_reduce(u64x4 x, u64x4 mod, u64x4 imod) {
#ifndef __AVX2__
        auto x_ninv = _mm256_mul_epu32(__m256i(x), __m256i(imod));
        auto x_res = _mm256_add_epi64(__m256i(x), _mm256_mul_epu32(x_ninv, __m256i(mod)));
        return u64x4(_mm256_bsrli_epi128(x_res, 4));
#else

        auto x_ninv = u64x4(u32x8(x) * u32x8(imod));
        return (x + x_ninv * mod) >> 32;
#endif
    }

    u64x4 montgomery_mul(u64x4 x, u64x4 y, u64x4 mod, u64x4 imod) {
#ifdef __AVX2__
        return montgomery_reduce(u64x4(_mm256_mul_epu32(__m256i(x), __m256i(y))), mod, imod);
#else
        return montgomery_reduce(x * y, mod, imod);
#endif
    }

    dx4 rotate_right(dx4 x) {
#ifdef __AVX2__
        return _mm256_permute4x64_pd(x, _MM_SHUFFLE(2, 1, 0, 3));
#else
        return __builtin_shufflevector(x, x, 3, 0, 1, 2);
#endif
    }
}
#endif // CP_ALGO_UTIL_SIMD_HPP
