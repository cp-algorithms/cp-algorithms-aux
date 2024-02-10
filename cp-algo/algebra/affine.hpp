#ifndef CP_ALGO_ALGEBRA_AFFINE_HPP
#define CP_ALGO_ALGEBRA_AFFINE_HPP
#include <optional>
#include <cassert>
namespace cp_algo::algebra {
    template<typename base>
    // a * x + b
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
            return lin(a * t.b + b * t.a, b * t.b + a * t.a * (*c), *c);
        }

        // a * (t.a * x + t.b) + b
        lin compose(lin const& t) const {
            return lin{a * t.a, a * t.b + b};
        }

        void prepend(lin const& t) {
            *this = t.compose(*this);
        }

        base eval(base x) const {
            return a * x + b;
        }
    };
}
#endif // CP_ALGO_ALGEBRA_AFFINE_HPP