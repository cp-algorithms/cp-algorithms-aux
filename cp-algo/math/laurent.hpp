#ifndef CP_ALGO_MATH_LAURENT_HPP
#define CP_ALGO_MATH_LAURENT_HPP
#include "fps.hpp"
namespace cp_algo::math {
    // x^offset times an FPS. offset is a lower bound, not a searched-for valuation.
    template<typename T>
    struct laurent {
        fps<T> series;
        int64_t offset = 0;
        T operator[](int64_t k) const {
            return k < offset ? T(0) : series[uint64_t(k) - uint64_t(offset)];
        }
        // Materialize n coefficients starting at the stored lower bound.
        poly_t<T> prefix(size_t n) const {return series.prefix(n);}
        friend laurent operator+(laurent a, laurent b) {
            int64_t lo = std::min(a.offset, b.offset);
            auto align = [lo](laurent p) {
                uint64_t shift = uint64_t(p.offset) - uint64_t(lo);
                return fps<T>([p, shift](size_t n, auto const&) {
                    return n < shift ? T(0) : p.series[n - shift];
                });
            };
            return {align(a) + align(b), lo};
        }
        friend laurent operator-(laurent a) {return {-a.series, a.offset};}
        friend laurent operator-(laurent a, laurent b) {return a + -b;}
        friend laurent operator*(laurent a, laurent b) {
            return {a.series * b.series, a.offset + b.offset};
        }
    };
    template<typename T>
    laurent<T> shift(laurent<T> p, int64_t k) {p.offset += k; return p;}

    // The coefficient at offset must be nonzero; no unbounded search for a leading term.
    template<typename T>
    laurent<T> inv(laurent<T> p) {return {inv(p.series), -p.offset};}

    template<typename T>
    laurent<T> operator/(laurent<T> a, laurent<T> b) {return a * inv(b);}

    template<typename T>
    laurent<T> deriv(laurent<T> p) {
        auto q = fps<T>([p](size_t n, auto const&) {return (T(p.offset) + T(n)) * p.series[n];});
        return {q, p.offset - 1};
    }
    // A Laurent antiderivative exists only when the x^-1 coefficient is zero.
    template<typename T>
    laurent<T> integr(laurent<T> p) {
        auto q = fps<T>([p](size_t n, auto const&) {
            T c = p.series[n], k = T(p.offset) + T(n) + T(1);
            if(k == T(0)) {
                if(c != T(0)) {throw std::domain_error("Laurent integral needs a logarithmic term");}
                return T(0);
            }
            return c / k;
        });
        return {q, p.offset + 1};
    }
}
#endif // CP_ALGO_MATH_LAURENT_HPP
