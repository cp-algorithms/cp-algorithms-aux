#ifndef CP_ALGO_ALGEBRA_AFFINE_HPP
#define CP_ALGO_ALGEBRA_AFFINE_HPP
#include <optional>
#include <utility>
#include <cassert>
namespace cp_algo::algebra {
    // a * x + b
    template<typename base>
    struct lin {
        base a = 1, b = 0;
        std::optional<base> c;
        lin() {}
        lin(base b): a(0), b(b) {}
        lin(base a, base b): a(a), b(b) {}
        lin(base a, base b, base _c): a(a), b(b), c(_c) {}

        // polynomial product modulo x^2 - c
        lin operator * (const lin& t) {
            assert(c && t.c && *c == *t.c);
            return {a * t.b + b * t.a, b * t.b + a * t.a * (*c), *c};
        }

        // a * (t.a * x + t.b) + b
        lin apply(lin const& t) const {
            return {a * t.a, a * t.b + b};
        }

        void prepend(lin const& t) {
            *this = t.apply(*this);
        }

        base eval(base x) const {
            return a * x + b;
        }
    };

    // (ax+b) / (cx+d)
    template<typename base>
    struct linfrac {
        base a, b, c, d;
        linfrac(): a(1), b(0), c(0), d(1) {} // x, identity for composition
        linfrac(base a): a(a), b(1), c(1), d(0) {} // a + 1/x, for continued fractions
        linfrac(base a, base b, base c, base d): a(a), b(b), c(c), d(d) {}
        
        // composition of two linfracs
        linfrac operator *(linfrac const& t) {
            auto [A, C] = apply(t.a, t.c);
            auto [B, D] = apply(t.b, t.d);
            return {A, B, C, D};
        }
        
        linfrac adj() {
            return {d, -b, -c, a};
        }
        
        // apply linfrac to A/B
        auto apply(base A, base B) {
            return std::pair{a * A + b * B, c * A + d * B};
        }
    };
}
#endif // CP_ALGO_ALGEBRA_AFFINE_HPP