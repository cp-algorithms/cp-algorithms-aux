#ifndef CP_ALGO_MATH_FPS_HPP
#define CP_ALGO_MATH_FPS_HPP
#include "poly/base.hpp"
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::fps_detail {
    // Tile pairs of positive indices by the highest power of two in their minimum.
    // A completed block contributes only to future coefficients: O(n log^2 n) in total.
    template<typename T>
    struct relaxed_product {
        big_vector<T> a, b, c;
        T pending(size_t n) {
            if(c.size() <= n) {c.resize(n + 1);}
            return c[n];
        }
        T append(T an, T bn) {
            size_t n = a.size();
            T res = pending(n) + an * (n ? b[0] : bn);
            if(n) {res += a[0] * bn;}
            a.push_back(an);
            b.push_back(bn);
            size_t end = n + 1;
            for(size_t len = 1; end % len == 0 && 2 * len <= end; len *= 2) {
                size_t start = end - len;
                c.resize(std::max(c.size(), end + 2 * len - 1));
                auto add = [&](auto const& x, auto const& y) {
                    if constexpr(modint_type<T>) {
                        if(len >= magic) {
                            big_vector<T> block(begin(x) + start, begin(x) + end);
                            fft::mul_truncate(block, std::span(y).subspan(len, len), 2 * len - 1);
                            for(size_t i = 0; i < block.size(); i++) {c[end + i] += block[i];}
                            return;
                        }
                    }
                    for(size_t i = 0; i < len; i++) {
                        for(size_t j = 0; j < len; j++) {c[end + i + j] += x[start + i] * y[len + j];}
                    }
                };
                add(a, b);
                if(start != len) {add(b, a);}
            }
            return res;
        }
    };
    template<typename T> struct product_transforms {};
    template<modint_type T> struct product_transforms<T> {
        std::vector<std::optional<fft::dft<T>>> transforms;
    };
    // One fixed input: each completed dyadic block fills the following block.
    // Cache fixed prefix transforms; only the middle of each product is needed.
    template<typename T, bool weighted = false>
    struct semi_relaxed_product: product_transforms<T> {
        std::span<T const> a;
        big_vector<T> b, c;
        T a0 = T(0);
        bool dirty = false;
        explicit semi_relaxed_product(std::span<T const> a): a(a) {}
        T coefficient(size_t n) const {
            return n < a.size() ? (weighted ? T(n) * a[n] : a[n]) : T(0);
        }
        // Exclude the current endpoint terms, matching relaxed_product::pending.
        T pending(size_t n) {
            flush();
            if(c.size() <= n) {c.resize(n + 1);}
            return c[n] - (n ? coefficient(n) * b[0] : T(0));
        }
        T append(T an, T bn) {
            size_t n = b.size();
            if(!n) {a0 = an;}
            T res = pending(n) + an * (n ? b[0] : bn);
            if(n) {res += a0 * bn;}
            b.push_back(bn);
            dirty = true;
            return res;
        }
        // Defer a completed block until its first output coefficient is requested.
        void flush() {
            if(!dirty) {return;}
            dirty = false;
            size_t end = b.size();
            if(a.size() < 2) {return;}
            size_t len = std::min(end & -end, std::bit_ceil(a.size() - 1));
            size_t start = end - len, count = std::min(2 * len, a.size());
            c.resize(std::max(c.size(), end + len));
            if constexpr(modint_type<T>) {
                if(len >= magic) {
                    size_t k = std::countr_zero(len);
                    if(this->transforms.size() <= k) {this->transforms.resize(k + 1);}
                    auto &fixed = this->transforms[k];
                    if(!fixed) {
                        big_vector<T> block(count);
                        for(size_t j = 0; j < count; j++) {block[j] = coefficient(j);}
                        fixed.emplace(block, len);
                    }
                    auto block = fft::dft<T>(std::span(b).subspan(start, len), len);
                    big_vector<T> result(2 * len);
                    // Retain a prefix only when reused; consume the first transform in place.
                    if(end == len) {
                        block.mul_inplace(*fixed, result, 2 * len);
                        fixed.reset();
                    } else {
                        block.mul(*fixed, result, 2 * len);
                    }
                    // Only the upper half is needed; cyclic wraparound affects the lower half.
                    for(size_t i = 0; i < len; i++) {c[end + i] += result[len + i];}
                    return;
                }
            }
            for(size_t i = 0; i < len; i++) {
                for(size_t j = len - i; j < std::min(count, 2 * len - i); j++) {
                    c[start + i + j] += b[start + i] * coefficient(j);
                }
            }
        }
    };
}
namespace cp_algo::math {
    // Shared, memoized coefficients. A generator may read its already computed prefix.
    // Copies share a cache; evaluation is single-threaded and must be causal.
    template<typename T>
    class fps {
        struct state {
            std::function<T(size_t, big_vector<T> const&)> generate;
            big_vector<T> cache;
            enum status {idle, fixed, active} mode = idle;
            T get(size_t n) {
                if(n < cache.size()) {return cache[n];}
                if(mode != idle) {
                    if(mode == fixed) {return T(0);}
                    throw std::logic_error("non-causal FPS dependency");
                }
                mode = active;
                try {
                    while(cache.size() <= n) {cache.push_back(generate(cache.size(), cache));}
                } catch(...) {mode = idle; throw;}
                mode = idle;
                return cache[n];
            }
        };
        std::shared_ptr<state> data;
        template<typename U> friend fps<U> inv(fps<U>);
        template<typename U> friend fps<U> log(fps<U>);
        template<typename U> friend fps<U> exp(fps<U>);
        template<bool weighted = false, typename F>
        fps with_product(F make) const {
            if(data->mode == state::fixed) {
                return make(fps_detail::semi_relaxed_product<T, weighted>(data->cache));
            }
            return make(fps_detail::relaxed_product<T>{});
        }
    public:
        fps(): fps(T(0)) {}
        fps(T c): fps(poly_t<T>(c)) {}
        fps(poly_t<T> p): data(std::make_shared<state>(
            [](size_t, auto const&) {return T(0);}, std::move(p.a), state::fixed)) {}
        template<typename F> requires std::invocable<F&, size_t, big_vector<T> const&>
        explicit fps(F f): data(std::make_shared<state>(std::move(f))) {}
        T operator[](size_t n) const {return data->get(n);}
        poly_t<T> prefix(size_t n) const {
            if(n) {data->get(n - 1);}
            big_vector<T> result(begin(data->cache), begin(data->cache) + std::min(n, data->cache.size()));
            result.resize(n);
            return result;
        }
        friend fps operator+(fps a, fps b) {
            return fps([a, b](size_t n, auto const&) {return a[n] + b[n];});
        }
        friend fps operator-(fps a, fps b) {
            return fps([a, b](size_t n, auto const&) {return a[n] - b[n];});
        }
        friend fps operator-(fps a) {
            return fps([a](size_t n, auto const&) {return -a[n];});
        }
        friend fps operator*(fps a, fps b) {
            if(a.data->mode != state::fixed && b.data->mode == state::fixed) {std::swap(a, b);}
            return a.with_product([a, b](auto product) {
                return fps([a, b, product = std::move(product)](size_t n, auto const&) mutable {
                    return product.append(a[n], b[n]);
                });
            });
        }
    };
    // Inverse with an invertible constant coefficient, evaluated on demand.
    template<typename T>
    fps<T> inv(fps<T> p) {
        return p.with_product([p](auto product) {
            return fps<T>([p, product = std::move(product), c = T(0)](size_t n, auto const&) mutable {
                if(n == 0) {
                    if(p[0] == T(0)) {throw std::domain_error("FPS inverse needs a nonzero constant");}
                    c = T(1) / p[0];
                    product.append(0, c);
                    return c;
                }
                T pn = p[n];
                T res = -(product.pending(n) + pn * c) * c;
                product.append(pn, res);
                return res;
            });
        });
    }
    template<typename T>
    fps<T> operator/(fps<T> a, fps<T> b) {return a * inv(b);}

    template<typename T>
    fps<T> deriv(fps<T> p) {
        return fps<T>([p](size_t n, auto const&) {return T(n + 1) * p[n + 1];});
    }
    template<typename T>
    fps<T> integr(fps<T> p) {
        return fps<T>([p](size_t n, auto const&) {return n ? p[n - 1] / T(n) : T(0);});
    }
    template<typename T>
    fps<T> log(fps<T> p) {
        // From p q' = p'; the product stores coefficient n of x q' as n*q[n].
        return p.with_product([p](auto product) {
            return fps<T>([p, product = std::move(product)](size_t n, auto const&) mutable {
                if(n == 0) {
                    if(p[0] != T(1)) {throw std::domain_error("FPS logarithm needs constant 1");}
                    product.append(0, 0);
                    return T(0);
                }
                T pn = p[n], qn = T(n) * pn - product.pending(n);
                product.append(pn, qn);
                return qn / T(n);
            });
        });
    }
    // From q' = p' q; completed convolution blocks are reused across requests.
    template<typename T>
    fps<T> exp(fps<T> p) {
        return p.template with_product<true>([p](auto product) {
            return fps<T>([p, product = std::move(product)](size_t n, auto const&) mutable {
                if(n == 0) {
                    if(p[0] != T(0)) {throw std::domain_error("FPS exponential needs constant 0");}
                    product.append(0, 1);
                    return T(1);
                }
                T pn = T(n) * p[n];
                T res = (product.pending(n) + pn) / T(n);
                product.append(pn, res);
                return res;
            });
        });
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_FPS_HPP
