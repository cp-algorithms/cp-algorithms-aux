#ifndef CP_ALGO_MATH_CVECTOR_HPP
#define CP_ALGO_MATH_CVECTOR_HPP
#include "../util/simd.hpp"
#include "../util/complex.hpp"
#include "../util/checkpoint.hpp"
#include "../util/big_alloc.hpp"
#include <immintrin.h>
#include <numbers>
#include <cstring>
#include <ranges>
#include <cmath>
#include <bit>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace stdx = std::experimental;
namespace cp_algo::math::fft {
    static constexpr size_t flen = 4;
    using ftype = double;
    using vftype = dx4;
    using point = complex<ftype>;
    using vpoint = complex<vftype>;
    static constexpr vftype vz = {};
    vpoint vi(vpoint const& r) {
        return {-imag(r), real(r)};
    }

    // Spectrum storage skips zero-filling on resize; every user writes before reading.
    template<class T>
    struct spectrum_alloc: big_alloc<T> {
        using big_alloc<T>::big_alloc;
        template<class U> struct rebind { using other = spectrum_alloc<U>; };
        template<class U> requires (std::is_same_v<U, vpoint>)
        void construct(U*) noexcept {}
    };
    using spectrum_vector = std::vector<vpoint, spectrum_alloc<vpoint>>;

    struct cvector {
        spectrum_vector r;
        cvector(size_t n) {
            n = std::max(flen, std::bit_ceil(n));
            r.assign(n / flen, vpoint{});
            prepare_roots(n / 16);
            checkpoint("cvector create");
        }

        vpoint& at(size_t k) {return r[k / flen];}
        vpoint at(size_t k) const {return r[k / flen];}
        template<class pt = point>
        inline void set(size_t k, pt const& t) {
            if constexpr(std::is_same_v<pt, point>) {
                real(r[k / flen])[k % flen] = real(t);
                imag(r[k / flen])[k % flen] = imag(t);
            } else {
                at(k) = t;
            }
        }
        template<class pt = point>
        inline pt get(size_t k) const {
            if constexpr(std::is_same_v<pt, point>) {
                return {real(r[k / flen])[k % flen], imag(r[k / flen])[k % flen]};
            } else {
                return at(k);
            }
        }

        size_t size() const {
            return flen * r.size();
        }
        static constexpr size_t eval_arg(size_t n) {
            if(n < pre_evals) {
                return eval_args[n];
            } else {
                return eval_arg(n / 2) | (n & 1) << (std::bit_width(n) - 1);
            }
        }
        static constexpr point eval_point(size_t n) {
            if(n % 2) {
                return -eval_point(n - 1);
            } else if(n % 4) {
                return eval_point(n - 2) * point(0, 1);
            } else if(n / 4 < pre_evals) {
                return evalp[n / 4];
            } else if(n / 4 - pre_evals < extra.size()) {
                return extra[n / 4 - pre_evals];
            } else {
                return polar<ftype>(1., std::numbers::pi / (ftype)std::bit_floor(n) * (ftype)eval_arg(n));
            }
        }
        static constexpr std::array<point, 32> roots = []() {
            std::array<point, 32> res;
            for(size_t i = 2; i < 32; i++) {
                res[i] = polar<ftype>(1., std::numbers::pi / (1ull << (i - 2)));
            }
            return res;
        }();
        static constexpr point root(size_t n) {
            return roots[std::bit_width(n)];
        }
        template<int step>
        static void exec_on_eval(size_t n, size_t k, auto &&callback) {
            callback(k, root(4 * step * n) * eval_point(step * k));
        }
        template<int step>
        static void exec_on_evals(size_t n, auto &&callback) {
            point factor=root(4*step*n);
            if constexpr(step==1 || step==2 || step==4){
                prepare_roots((step*n+3)/4);
                size_t i=0;
                if constexpr(step==1){
                    for(;i+4<=n;i+=4){
                        size_t k=i/4;
                        point e=k<pre_evals?evalp[k]:extra[k-pre_evals];
                        point v=factor*e;
                        callback(i,v);callback(i+1,-v);
                        point iv(-imag(v),real(v));
                        callback(i+2,iv);callback(i+3,-iv);
                    }
                }else if constexpr(step==2){
                    for(;i+2<=n;i+=2){
                        size_t k=i/2;
                        point e=k<pre_evals?evalp[k]:extra[k-pre_evals];
                        point v=factor*e;
                        callback(i,v);callback(i+1,point(-imag(v),real(v)));
                    }
                }
                for(;i<n;i++){
                    size_t index=step*i,k=index/4;
                    point e=k<pre_evals?evalp[k]:extra[k-pre_evals];
                    if(index&2)e=e*point(0,1);
                    if(index&1)e=-e;
                    callback(i,factor*e);
                }
            }else{
                for(size_t i=0;i<n;i++)callback(i,factor*eval_point(step*i));
            }
        }

        static void do_dot_iter(point rt, vpoint& Bv, vpoint const& Av, vpoint& res) {
            res += Av * Bv;
            real(Bv) = rotate_right(real(Bv));
            imag(Bv) = rotate_right(imag(Bv));
            auto x = real(Bv)[0], y = imag(Bv)[0];
            real(Bv)[0] = x * real(rt) - y * imag(rt);
            imag(Bv)[0] = x * imag(rt) + y * real(rt);
        }

        template<size_t fixed = 0>
        void dot(cvector const& t) {
            size_t n = fixed?fixed:this->size();
            exec_on_evals<1>(n / flen, [&](size_t k, point rt) __attribute__((always_inline)) {
                k *= flen;
                auto [Ax, Ay] = at(k);
                auto Cv = t.at(k);
                vpoint vrt = {vz + real(rt), vz + imag(rt)};
                auto Cr = Cv * vrt;
                vpoint res = vz;
                auto iter = [&]<int i>() __attribute__((always_inline)) {
                    auto wrap = [&](vftype original, vftype rotated) {
                        if constexpr(i == 0) {return original;}
                        else {return __builtin_shufflevector(rotated, original, 4-i,5-i,6-i,7-i);}
                    };
                    vpoint Cw = {wrap(real(Cv),real(Cr)),wrap(imag(Cv),imag(Cr))};
                    vpoint Av = {vz+Ax[i],vz+Ay[i]};
                    return Av*Cw;
                };
                auto p0=iter.template operator()<0>(), p1=iter.template operator()<1>();
                auto p2=iter.template operator()<2>(), p3=iter.template operator()<3>();
                res=(p0+p1)+(p2+p3);
                set(k, res);
            });
            checkpoint("dot");
        }
        // normalize=false leaves the inverse-transform scale for the caller.
        template<bool partial = true, bool normalize = true, size_t fixed = 0>
        void ifft() {
            size_t n = fixed?fixed:size();
            if constexpr (!partial) {
                prepare_roots(n / 4);
                point pi(0, 1);
                exec_on_evals<4>(n / 4, [&](size_t k, point rt) __attribute__((always_inline)) {
                    k *= 4;
                    point v1 = conj(rt);
                    point v2 = v1 * v1;
                    point v3 = v1 * v2;
                    auto A = get(k);
                    auto B = get(k + 1);
                    auto C = get(k + 2);
                    auto D = get(k + 3);
                    set(k, (A + B) + (C + D));
                    set(k + 2, ((A + B) - (C + D)) * v2);
                    set(k + 1, ((A - B) - pi * (C - D)) * v1);
                    set(k + 3, ((A - B) + pi * (C - D)) * v3);
                });
            }
            bool parity = std::countr_zero(n) % 2;
            if(parity) {
                exec_on_evals<2>(n / (2 * flen), [&](size_t k, point rt) __attribute__((always_inline)) {
                    k *= 2 * flen;
                    vpoint cvrt = {vz + real(rt), vz - imag(rt)};
                    auto B = at(k) - at(k + flen);
                    at(k) += at(k + flen);
                    at(k + flen) = B * cvrt;
                });
            }

            transform<true,fixed>(n, parity);
            checkpoint("ifft");
            if constexpr(normalize) {
                auto scale = vz + ftype(partial ? flen : 1) / ftype(n);
                for(size_t k = 0; k < n; k += flen) {
                    set(k, get<vpoint>(k) * scale);
                }
            }
        }
        template<bool partial = true, size_t fixed = 0>
        void fft() {
            size_t n = fixed?fixed:size();
            prepare_roots(n / (partial ? 16 : 4));
            bool parity = std::countr_zero(n) % 2;
            transform<false,fixed>(n, parity);
            if(parity) {
                exec_on_evals<2>(n / (2 * flen), [&](size_t k, point rt) __attribute__((always_inline)) {
                    k *= 2 * flen;
                    vpoint vrt = {vz + real(rt), vz + imag(rt)};
                    auto t = at(k + flen) * vrt;
                    at(k + flen) = at(k) - t;
                    at(k) += t;
                });
            }
            if constexpr (!partial) {
                prepare_roots(n / 4);
                point pi(0, 1);
                exec_on_evals<4>(n / 4, [&](size_t k, point rt) __attribute__((always_inline)) {
                    k *= 4;
                    point v1 = rt;
                    point v2 = v1 * v1;
                    point v3 = v1 * v2;
                    auto A = get(k);
                    auto B = get(k + 1) * v1;
                    auto C = get(k + 2) * v2;
                    auto D = get(k + 3) * v3;
                    set(k, (A + C) + (B + D));
                    set(k + 1, (A + C) - (B + D));
                    set(k + 2, (A - C) + pi * (B - D));
                    set(k + 3, (A - C) - pi * (B - D));
                });
            }
            checkpoint("fft");
        }
        static std::array<vpoint,4> transpose(std::array<vpoint,4> const&a){
            auto half=[](auto part){
                auto a0=__m256d(part(0)),a1=__m256d(part(1)),a2=__m256d(part(2)),a3=__m256d(part(3));
                auto t0=_mm256_unpacklo_pd(a0,a1),t1=_mm256_unpackhi_pd(a0,a1),t2=_mm256_unpacklo_pd(a2,a3),t3=_mm256_unpackhi_pd(a2,a3);
                return std::array<vftype,4>{vftype(_mm256_permute2f128_pd(t0,t2,0x20)),vftype(_mm256_permute2f128_pd(t1,t3,0x20)),vftype(_mm256_permute2f128_pd(t0,t2,0x31)),vftype(_mm256_permute2f128_pd(t1,t3,0x31))};
            };
            auto re=half([&](int k){return real(a[k]);}),im=half([&](int k){return imag(a[k]);});
            return {vpoint{re[0],im[0]},vpoint{re[1],im[1]},vpoint{re[2],im[2]},vpoint{re[3],im[3]}};
        }
        template<size_t fixed>void dot_fused16(cvector const& t,size_t offset,size_t length){
            constexpr size_t n=fixed;
            point factor=root(n);
            for(size_t pos=offset;pos<offset+length;pos+=16){
                size_t k=pos/16;point e=k<pre_evals?evalp[k]:extra[k-pre_evals];point rt=factor*e;
                vpoint v1={vz+real(rt),vz+imag(rt)},v2=v1*v1,v3=v1*v2;
                auto forward=[&](cvector const& a) __attribute__((always_inline)) {
                    auto A=a.at(pos),B=a.at(pos+4)*v1,C=a.at(pos+8)*v2,D=a.at(pos+12)*v3;
                    return std::array<vpoint,4>{(A+C)+(B+D),(A+C)-(B+D),(A-C)+vi(B-D),(A-C)-vi(B-D)};
                };
                auto a=transpose(forward(*this)),b=transpose(forward(t));
                vpoint w={vftype{real(rt),-real(rt),-imag(rt),imag(rt)},vftype{imag(rt),-imag(rt),real(rt),-real(rt)}};
                auto cmadd=[](vpoint x,vpoint y,vpoint c) __attribute__((always_inline)) {
                    auto re=_mm256_fmadd_pd(__m256d(real(x)),__m256d(real(y)),_mm256_fnmadd_pd(__m256d(imag(x)),__m256d(imag(y)),__m256d(real(c))));
                    auto im=_mm256_fmadd_pd(__m256d(real(x)),__m256d(imag(y)),_mm256_fmadd_pd(__m256d(imag(x)),__m256d(real(y)),__m256d(imag(c))));
                    return vpoint{vftype(re),vftype(im)};
                };
                auto cmcross=[](vpoint x,vpoint y,vpoint p,vpoint q) __attribute__((always_inline)) {
                    auto re=_mm256_sub_pd(_mm256_fmsub_pd(__m256d(real(x)),__m256d(real(y)),__m256d(real(p))),_mm256_fmadd_pd(__m256d(imag(x)),__m256d(imag(y)),__m256d(real(q))));
                    auto im=_mm256_add_pd(_mm256_fmsub_pd(__m256d(real(x)),__m256d(imag(y)),__m256d(imag(p))),_mm256_fmsub_pd(__m256d(imag(x)),__m256d(real(y)),__m256d(imag(q))));
                    return vpoint{vftype(re),vftype(im)};
                };
                auto mul=[&](vpoint a0,vpoint a1,vpoint b0,vpoint b1) __attribute__((always_inline)) {
                    auto p=a0*b0,q=a1*b1;
                    return std::array<vpoint,2>{cmadd(w,q,p),cmcross(a0+a1,b0+b1,p,q)};
                };
                auto p=mul(a[0],a[2],b[0],b[2]),q=mul(a[1],a[3],b[1],b[3]);
                auto m=mul(a[0]+a[1],a[2]+a[3],b[0]+b[1],b[2]+b[3]);
                auto c=transpose({cmadd(w,q[1],p[0]),(m[0]-p[0])-q[0],p[1]+q[0],(m[1]-p[1])-q[1]});
                auto A=c[0],B=c[1],C=c[2],D=c[3];
                // Recompute the inverse weights in the original operation order.
                vpoint u1={vz+real(rt),vz-imag(rt)},u2=u1*u1,u3=u1*u2;
                at(pos)=(A+B)+(C+D);at(pos+8)=((A+B)-(C+D))*u2;
                at(pos+4)=((A-B)-vi(C-D))*u1;at(pos+12)=((A-B)+vi(C-D))*u3;
            }
        }
        // Radix-64 out-of-cache pass for n = 2^24, run as two tiled radix-8 stages.
        // With Fuse=1 (forward only) the first stage lifts u32 residues to Gaussian
        // coordinates on load, so the spectrum buffer is written once and never re-read
        // before the in-cache block phase. The index-keyed noise makes both branches see
        // the same representatives.
        struct fuse_args {
            const uint32_t* src = nullptr;
            size_t count = 0;
            uint64_t seed = 0;
            double a = 0, b = 0, a_over_p = 0, b_over_p = 0;
        };
        static constexpr size_t sweep_tile = 256;
        template<bool inverse, int Mode, size_t Tile, int Fuse = 0, bool Neg = false>
        void sweep8(fuse_args const& fa = {}) {
            constexpr size_t n = 1 << 24;
            static const std::array<std::array<point, 7>, 9> weights = []() {
                std::array<std::array<point, 7>, 9> table;
                for(size_t count: {size_t(1), size_t(8)}) for(size_t k = 0; k < count; k++) {
                    size_t rev = 0, x = k;
                    for(size_t c = count; c > 1; c >>= 1) {rev = (rev << 1) | (x & 1); x >>= 1;}
                    long double angle = std::numbers::pi_v<long double> * (1 + 4 * rev) / (16 * count);
                    if(k & 1) {angle -= std::numbers::pi_v<long double> / 4;}
                    for(size_t j = 1; j < 8; j++) {
                        long double a = j * angle;
                        if constexpr(Mode == 0) {table[(count - 1) / 7 + k][j - 1] = {double(cosl(a)), double(sinl(a))};}
                        else {table[(count - 1) / 7 + k][j - 1] = {double(sinl(a)), double(tanl(a / 2))};}
                    }
                }
                return table;
            }();
            auto lift = [&](size_t idx) __attribute__((always_inline)) -> vpoint {
                i32x4 bits{};
                if(idx + 4 <= fa.count) {std::memcpy(&bits, fa.src + idx, sizeof(bits));}
                else if(idx < fa.count) {for(size_t j = 0; j < fa.count - idx; j++) {bits[j] = int32_t(fa.src[idx + j]);}}
                else {return vpoint{vz, vz};}
                auto x = __builtin_convertvector(bits, vftype);
                u32x4 h = u32x4{uint32_t(idx), uint32_t(idx + 1), uint32_t(idx + 2), uint32_t(idx + 3)} ^ uint32_t(fa.seed);
                h *= 0x9E3779B1u; h ^= h >> 15; h *= 0x85EBCA77u; h ^= h >> 13; h *= 0xC2B2AE3Du; h ^= h >> 16;
                auto noise = __builtin_convertvector(i32x4(h), vftype) * 0x1p-32;
                auto q = round(x * fa.a_over_p + noise), t = round(x * fa.b_over_p + noise);
                auto re = x - q * fa.a - t * fa.b, im = t * fa.a - q * fa.b;
                return vpoint{re, Neg ? -im : im};
            };
            auto stage = [&]<bool top>(size_t offset, size_t length, size_t begin, size_t end) __attribute__((always_inline)) {
                size_t step = length / 8, k = offset / length;
                auto const& w = weights[(n / length - 1) / 7 + k];
                auto rot = [&]<size_t J>(vpoint z) __attribute__((always_inline)) {
                    auto c = w[J - 1];
                    if constexpr(Mode == 0) {return z * vpoint{vz + real(c), inverse ? vz - imag(c) : vz + imag(c)};}
                    else {
                        auto s = inverse ? vz - real(c) : vz + real(c), t = inverse ? vz - imag(c) : vz + imag(c);
                        auto x = vftype(_mm256_fnmadd_pd(__m256d(t), __m256d(imag(z)), __m256d(real(z))));
                        auto y = vftype(_mm256_fmadd_pd(__m256d(s), __m256d(x), __m256d(imag(z))));
                        return vpoint{vftype(_mm256_fnmadd_pd(__m256d(t), __m256d(y), __m256d(x))), y};
                    }
                };
                auto add = [](vpoint a, vpoint b) __attribute__((always_inline)) {
                    if constexpr(Mode != 2) {return a + b;}
                    else {return vpoint{vftype(_mm256_fmadd_pd(__m256d(real(a)), _mm256_set1_pd(1), __m256d(real(b)))), vftype(_mm256_fmadd_pd(__m256d(imag(a)), _mm256_set1_pd(1), __m256d(imag(b))))};}
                };
                auto sub = [](vpoint a, vpoint b) __attribute__((always_inline)) {
                    if constexpr(Mode != 2) {return a - b;}
                    else {return vpoint{vftype(_mm256_fmsub_pd(__m256d(real(a)), _mm256_set1_pd(1), __m256d(real(b)))), vftype(_mm256_fmsub_pd(__m256d(imag(a)), _mm256_set1_pd(1), __m256d(imag(b))))};}
                };
                auto d4 = [](vpoint a, vpoint b, vpoint c, vpoint d) __attribute__((always_inline)) {
                    auto s = a + c, t = a - c, u = b + d, v = vi(b - d);
                    if constexpr(inverse) {return std::array<vpoint, 4>{s + u, t - v, s - u, t + v};}
                    else {return std::array<vpoint, 4>{s + u, t + v, s - u, t - v};}
                };
                constexpr double q = 0.707106781186547524400844362104849039;
                auto r1 = [](vpoint z) __attribute__((always_inline)) {
                    if constexpr(inverse) {return vpoint{(real(z) + imag(z)) * q, (imag(z) - real(z)) * q};}
                    else {return vpoint{(real(z) - imag(z)) * q, (real(z) + imag(z)) * q};}
                };
                auto r3 = [](vpoint z) __attribute__((always_inline)) {
                    if constexpr(inverse) {return vpoint{(imag(z) - real(z)) * q, (-imag(z) - real(z)) * q};}
                    else {return vpoint{(-real(z) - imag(z)) * q, (real(z) - imag(z)) * q};}
                };
                std::array<vpoint*, 8> input, output;
                for(size_t j = 0; j < 8; j++) {input[j] = output[j] = r.data() + (offset + begin + j * step) / flen;}
                if(k & 1) {
                    constexpr std::array<size_t, 8> perm = {7, 6, 4, 5, 0, 1, 2, 3};
                    for(size_t j = 0; j < 8; j++) {
                        if constexpr(inverse) {input[j] = output[perm[j]];}
                        else {output[j] = input[perm[j]];}
                    }
                }
                constexpr bool fused_in = top && Fuse == 1 && !inverse;
                for(size_t j = 0; j < (end - begin) / flen; j++) {
                    size_t base = offset + begin + j * flen;
                    auto in = [&]<size_t S>() __attribute__((always_inline)) {
                        if constexpr(fused_in) {return lift(base + S * step);}
                        else {return input[S][j];}
                    };
                    if constexpr(!inverse) {
                        auto E = d4(in.template operator()<0>(), rot.template operator()<2>(in.template operator()<2>()), rot.template operator()<4>(in.template operator()<4>()), rot.template operator()<6>(in.template operator()<6>()));
                        auto O = d4(rot.template operator()<1>(in.template operator()<1>()), rot.template operator()<3>(in.template operator()<3>()), rot.template operator()<5>(in.template operator()<5>()), rot.template operator()<7>(in.template operator()<7>()));
                        auto a = O[0], b = r1(O[1]), c = vi(O[2]), d = r3(O[3]);
                        output[0][j] = add(E[0], a); output[1][j] = sub(E[0], a);
                        output[2][j] = add(E[2], c); output[3][j] = sub(E[2], c);
                        output[4][j] = add(E[1], b); output[5][j] = sub(E[1], b);
                        output[6][j] = add(E[3], d); output[7][j] = sub(E[3], d);
                    } else {
                        auto E = d4(add(input[0][j], input[1][j]), add(input[4][j], input[5][j]), add(input[2][j], input[3][j]), add(input[6][j], input[7][j]));
                        auto O = d4(sub(input[0][j], input[1][j]), r1(sub(input[4][j], input[5][j])), -vi(sub(input[2][j], input[3][j])), r3(sub(input[6][j], input[7][j])));
                        output[0][j] = E[0]; output[1][j] = rot.template operator()<1>(O[0]);
                        output[2][j] = rot.template operator()<2>(E[1]); output[3][j] = rot.template operator()<3>(O[1]);
                        output[4][j] = rot.template operator()<4>(E[2]); output[5][j] = rot.template operator()<5>(O[2]);
                        output[6][j] = rot.template operator()<6>(E[3]); output[7][j] = rot.template operator()<7>(O[3]);
                    }
                }
            };
            constexpr size_t h = n / 64;
            for(size_t j = 0; j < h; j += Tile) {
                size_t end = std::min(h, j + Tile);
                auto first = [&]() {for(size_t t = 0; t < 8; t++) {stage.template operator()<true>(0, n, j + t * h, end + t * h);}};
                auto second = [&]() {for(size_t k = 0; k < 8; k++) {stage.template operator()<false>(k * n / 8, n / 8, j, end);}};
                if constexpr(inverse) {second(); first();} else {first(); second();}
            }
        }
        // Product for n = 2^24: the top three stages split each input into 64 independent
        // blocks. Both forward transforms read the u32 inputs directly; each block then
        // completes its two forward transforms, the product, and the inverse before the
        // final three inverse stages combine the results.
        template<bool Neg>
        void cache_product(cvector& b, fuse_args const& fa, fuse_args const& fb) {
            constexpr size_t n = 1 << 24, block = 1 << 18;
            prepare_roots(n / 16); prepare_shear_roots();
            sweep8<false, 2, sweep_tile, 1, Neg>(fa);
            b.sweep8<false, 2, sweep_tile, 1, Neg>(fb);
            checkpoint("fused forward");
            for(size_t offset = 0; offset < n; offset += block) {
                transform<false, n, block, 0, true>(n, false, offset, block);
                b.transform<false, n, block, 0, true>(n, false, offset, block);
                dot_fused16<n>(b, offset, block);
                transform<true, n, block, 0, true>(n, false, offset, block);
            }
            checkpoint("blocks");
            sweep8<true, 2, sweep_tile>();
            checkpoint("sweep inverse");
        }
        static constexpr size_t pre_evals = 1 << 16;
        static const std::array<size_t, pre_evals> eval_args;
        static const std::array<point, pre_evals> evalp;
    private:
        // Tile two radix-four stages together before descending into each child.
        template<bool inverse, size_t fixed = 0, size_t range_fixed=0, int top_fixed=-1,bool omit16=false>
        void transform(size_t input_n, bool parity, size_t range_offset=0, size_t range_length=0, int top_only=0) {
            if constexpr(range_fixed)range_length=range_fixed;
            if constexpr(top_fixed>=0)top_only=top_fixed;
            const size_t n=fixed?fixed:input_n;
            if constexpr(!range_fixed){prepare_roots(n/16);
            if constexpr(fixed==(1<<24))prepare_shear_roots();}
            size_t log_n=std::countr_zero(n);
            auto butterfly = [&](size_t offset,size_t length,size_t begin,size_t end) __attribute__((always_inline)) {
                if constexpr(omit16)if(length==16)return;
                size_t step=length/4,log_length=std::countr_zero(length),k=offset>>log_length;
                auto *p0=r.data()+(offset+begin)/flen,*p1=r.data()+(offset+begin+step)/flen,*p2=r.data()+(offset+begin+2*step)/flen,*p3=r.data()+(offset+begin+3*step)/flen;
                auto run=[&]<bool shear>() __attribute__((always_inline)) {
                    vpoint v1,v2,v3;vftype t1{},t2{},t3{};
                    if constexpr(shear && fixed==(1<<24)) {
                        auto const& c=shear_roots[((n>>log_length)-1)/3+k];
                        v1={vz,inverse?vz-real(c[0]):vz+real(c[0])};
                        v2={vz,inverse?vz-real(c[1]):vz+real(c[1])};
                        v3={vz,inverse?vz-real(c[2]):vz+real(c[2])};
                        t1=inverse?vz-imag(c[0]):vz+imag(c[0]);
                        t2=inverse?vz-imag(c[1]):vz+imag(c[1]);
                        t3=inverse?vz-imag(c[2]):vz+imag(c[2]);
                    }else {
                        point e=k<pre_evals?evalp[k]:extra[k-pre_evals];point rt=roots[log_n+5-log_length]*e;
                        v1={vz+real(rt),inverse?vz-imag(rt):vz+imag(rt)};
                        if constexpr(shear)if(k&1){if constexpr(inverse)v1=vi(v1);else v1=-vi(v1);}
                        v2=v1*v1;v3=v1*v2;
                        if constexpr(shear){t1=imag(v1)/(vz+1.0+real(v1));t2=imag(v2)/(vz+1.0+real(v2));t3=imag(v3)/(vz+1.0+real(v3));}
                    }
                    auto rotate=[](vpoint z,vpoint v,vftype t) __attribute__((always_inline)) {
                        if constexpr(!shear)return z*v;
                        else {
                            auto x=vftype(_mm256_fnmadd_pd(__m256d(t),__m256d(imag(z)),__m256d(real(z))));
                            auto y=vftype(_mm256_fmadd_pd(__m256d(imag(v)),__m256d(x),__m256d(imag(z))));
                            return vpoint{vftype(_mm256_fnmadd_pd(__m256d(t),__m256d(y),__m256d(x))),y};
                        }
                    };
                    auto *i0=p0,*i1=p1,*i2=p2,*i3=p3,*o0=p0,*o1=p1,*o2=p2,*o3=p3;
                    if constexpr(shear)if(k&1) {
                        if constexpr(inverse){i0=p3;i1=p2;i2=p0;i3=p1;}
                        else{o0=p3;o1=p2;o2=p0;o3=p1;}
                    }
                    for(size_t j=0;j<(end-begin)/flen;j++) {
                        auto A=i0[j],B=i1[j],C=i2[j],D=i3[j];
                        if constexpr(inverse) {
                            o0[j]=(A+B)+(C+D);
                            o2[j]=rotate((A+B)-(C+D),v2,t2);
                            o1[j]=rotate((A-B)-vi(C-D),v1,t1);
                            o3[j]=rotate((A-B)+vi(C-D),v3,t3);
                        }else{
                            B=rotate(B,v1,t1);C=rotate(C,v2,t2);D=rotate(D,v3,t3);
                            o0[j]=(A+C)+(B+D);o1[j]=(A+C)-(B+D);
                            o2[j]=(A-C)+vi(B-D);o3[j]=(A-C)-vi(B-D);
                        }
                    }
                };
                if(length>=256)run.template operator()<true>();else run.template operator()<false>();
            };
            if(top_only){
                size_t offset=range_offset,length=range_length;
                if(top_only==1){butterfly(offset,length,0,length/4);return;}
                if(top_only==3){
                    size_t h=length/64;
                    for(size_t j=0;j<h;j+=512){
                        size_t end=std::min(h,j+512);
                        auto stage0=[&](){for(size_t t=0;t<16;t++)butterfly(offset,length,j+t*h,end+t*h);};
                        auto stage1=[&](){for(size_t q=0;q<4;q++)for(size_t t=0;t<4;t++)butterfly(offset+q*length/4,length/4,j+t*h,end+t*h);};
                        auto stage2=[&](){for(size_t q=0;q<16;q++)butterfly(offset+q*length/16,length/16,j,end);};
                        if constexpr(inverse){stage2();stage1();stage0();}else{stage0();stage1();stage2();}
                    }
                    return;
                }
                size_t step=length/16;
                for(size_t j=0;j<step;j+=256){
                    size_t end=std::min(step,j+256);
                    if constexpr(inverse){
                        for(size_t t=0;t<4;t++)butterfly(offset+t*length/4,length/4,j,end);
                        for(size_t t=0;t<4;t++)butterfly(offset,length,j+t*step,end+t*step);
                    }else{
                        for(size_t t=0;t<4;t++)butterfly(offset,length,j+t*step,end+t*step);
                        for(size_t t=0;t<4;t++)butterfly(offset+t*length/4,length/4,j,end);
                    }
                }
                return;
            }
            auto recurse = [&](auto &&self, size_t offset, size_t length) -> void {
                if(length < 4 * flen) {return;}
                if(length >= (1 << 15)) {
                    size_t step = length / 16;
                    if constexpr(inverse) {
                        for(size_t t = 0; t < 16; t++) {self(self, offset + t*step, step);}
                    }
                    for(size_t j = 0; j < step; j += 256) {
                        size_t end = std::min(step, j+256);
                        if constexpr(inverse) {
                            for(size_t t=0;t<4;t++) {butterfly(offset+t*length/4, length/4, j,end);}
                            for(size_t t=0;t<4;t++) {butterfly(offset,length,j+t*step,end+t*step);}
                        } else {
                            for(size_t t=0;t<4;t++) {butterfly(offset,length,j+t*step,end+t*step);}
                            for(size_t t=0;t<4;t++) {butterfly(offset+t*length/4,length/4,j,end);}
                        }
                    }
                    if constexpr(!inverse) {
                        for(size_t t = 0; t < 16; t++) {self(self, offset + t*step, step);}
                    }
                } else if(length >= (size_t(1) << (6 + parity))) {
                    auto finish=[&]<bool par>() {
                        constexpr size_t chunk=size_t(1)<<(6+par);
                        constexpr size_t bottom=size_t(1)<<(4+par);
                        auto small=[&]<size_t L>(auto&& self,size_t pos) __attribute__((always_inline)) -> void {
                            if constexpr(inverse && L>bottom) {
                                for(size_t q=0;q<4;q++)self.template operator()<L/4>(self,pos+q*(L/4));
                            }
                            butterfly(pos,L,0,L/4);
                            if constexpr(!inverse && L>bottom) {
                                for(size_t q=0;q<4;q++)self.template operator()<L/4>(self,pos+q*(L/4));
                            }
                        };
                        for(size_t leaf=offset;leaf<offset+length;leaf+=chunk){
                            if constexpr(inverse){
                                small.template operator()<chunk>(small,leaf);
                                size_t level=std::min<size_t>(std::countr_one(leaf+chunk-1),std::countr_zero(length));
                                for(size_t lvl=6+2+par;lvl<=level;lvl+=2){size_t len=size_t(1)<<lvl;butterfly(leaf & ~(len-1),len,0,len/4);}
                            }else{
                                size_t level=std::min<size_t>(std::countr_zero(n+leaf),std::countr_zero(length));
                                level-=level%2!=par;
                                for(size_t lvl=level;lvl>=6+2+par;lvl-=2){size_t len=size_t(1)<<lvl;butterfly(leaf & ~(len-1),len,0,len/4);}
                                small.template operator()<chunk>(small,leaf);
                            }
                        }
                    };
                    if(parity)finish.template operator()<true>();else finish.template operator()<false>();
                } else {
                    if constexpr(inverse) {
                        for(size_t leaf = offset + 3 * flen; leaf < offset + length; leaf += 4 * flen) {
                            size_t level = std::min<size_t>(std::countr_one(leaf + 3), std::countr_zero(length));
                            for(size_t lvl = 4 + parity; lvl <= level; lvl += 2) {
                                size_t len = size_t(1) << lvl;
                                butterfly(leaf & ~(len-1), len, 0, len / 4);
                            }
                        }
                    } else {
                        for(size_t leaf = offset; leaf < offset + length; leaf += 4 * flen) {
                            size_t level = std::min<size_t>(std::countr_zero(n + leaf), std::countr_zero(length));
                            level -= level % 2 != parity;
                            for(size_t lvl = level; lvl >= 4; lvl -= 2) {
                                size_t len = size_t(1) << lvl;
                                butterfly(leaf & ~(len-1), len, 0, len / 4);
                            }
                        }
                    }
                }
            };
            // Radix two is performed separately at the leaves.
            recurse(recurse, range_offset, range_length?range_length:n);
        }
        static big_vector<std::array<point,3>> shear_roots;
        static void prepare_shear_roots() {
            constexpr size_t n=1<<24;
            if(!shear_roots.empty())return;
            shear_roots.resize((n/64-1)/3,std::array<point,3>{});
            for(size_t len=n;len>=256;len/=4) {
                size_t count=n/len,base=(count-1)/3;
                point factor=roots[29-std::countr_zero(len)];
                for(size_t k=0;k<count;k++) {
                    point rt=factor*(k<pre_evals?evalp[k]:extra[k-pre_evals]);
                    vpoint v1={vz+real(rt),vz+imag(rt)};
                    if(k&1)v1=-vi(v1);
                    vpoint v2=v1*v1,v3=v1*v2;
                    vftype t1=imag(v1)/(vz+1.0+real(v1)),t2=imag(v2)/(vz+1.0+real(v2)),t3=imag(v3)/(vz+1.0+real(v3));
                    shear_roots[base+k]={point(imag(v1)[0],t1[0]),point(imag(v2)[0],t2[0]),point(imag(v3)[0],t3[0])};
                }
            }
        }
        static big_vector<point> extra;
        // Keep the usual table small; cache additional roots for large transforms.
        static void prepare_roots(size_t n) {
            if(n <= pre_evals + extra.size()) {return;}
            size_t old = extra.size();
            extra.resize(std::bit_ceil(n) - pre_evals);
            static const std::array<point,256> coarse=[](){
                std::array<point,256> out;
                for(size_t i=0;i<256;i++)out[i]=polar<ftype>(1.,std::numbers::pi*double((eval_args[256+i]-1)/2)/512.0);
                return out;
            }();
            for(size_t h=pre_evals+old;h<pre_evals+extra.size();h*=2){
                for(size_t i=h;i<2*h;i+=256){
                    point fine=polar<ftype>(1.,std::numbers::pi/double(4*h)*double(eval_arg(4*i)));
                    for(size_t j=0;j<256;j++)extra[i+j-pre_evals]=coarse[j]*fine;
                }
            }
        }
    };

    big_vector<std::array<point,3>> cvector::shear_roots;
    big_vector<point> cvector::extra;

    const std::array<size_t, cvector::pre_evals> cvector::eval_args = []() {
        std::array<size_t, pre_evals> res = {};
        for(size_t i = 1; i < pre_evals; i++) {
            res[i] = res[i >> 1] | (i & 1) << (std::bit_width(i) - 1);
        }
        return res;
    }();
    const std::array<point, cvector::pre_evals> cvector::evalp = []() {
        std::array<point, pre_evals> res = {};
        res[0] = 1;
        for(size_t n = 1; n < pre_evals; n++) {
            res[n] = polar<ftype>(1., std::numbers::pi * ftype(eval_args[n]) / ftype(4 * std::bit_floor(n)));
        }
        return res;
    }();
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_CVECTOR_HPP
