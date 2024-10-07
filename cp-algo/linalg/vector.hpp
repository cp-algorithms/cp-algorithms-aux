#ifndef CP_ALGO_LINALG_VECTOR_HPP
#define CP_ALGO_LINALG_VECTOR_HPP
#include "../random/rng.hpp"
#include "../number_theory/modint.hpp"
#include <functional>
#include <algorithm>
#include <valarray>
#include <iostream>
#include <iterator>
#include <cassert>
namespace cp_algo::linalg {
    template<class vec, typename base>
    struct valarray_base: std::valarray<base> {
        using Base = std::valarray<base>;
        using Base::Base;

        valarray_base(base const& t): Base(t, 1) {}

        auto begin() {return std::begin(to_valarray());}
        auto begin() const {return std::begin(to_valarray());}
        auto end() {return std::end(to_valarray());}
        auto end() const {return std::end(to_valarray());}

        bool operator == (vec const& t) const {return std::ranges::equal(*this, t);}
        bool operator != (vec const& t) const {return !(*this == t);}

        vec operator-() const {return Base::operator-();}

        static vec from_range(auto const& R) {
            vec res(std::ranges::distance(R));
            std::ranges::copy(R, res.begin());
            return res;
        }
        Base& to_valarray() {return static_cast<Base&>(*this);}
        Base const& to_valarray() const {return static_cast<Base const&>(*this);}
    };

    template<class vec, typename base>
    vec operator+(valarray_base<vec, base> const& a, valarray_base<vec, base> const& b) {
        return a.to_valarray() + b.to_valarray();
    }
    template<class vec, typename base>
    vec operator-(valarray_base<vec, base> const& a, valarray_base<vec, base> const& b) {
        return a.to_valarray() - b.to_valarray();
    }

    template<class vec, typename base>
    struct vec_base: valarray_base<vec, base> {
        using Base = valarray_base<vec, base>;
        using Base::Base;

        static vec ei(size_t n, size_t i) {
            vec res(n);
            res[i] = 1;
            return res;
        }

        virtual void add_scaled(vec const& b, base scale, size_t i = 0) {
            if(scale != base(0)) {
                for(; i < size(*this); i++) {
                    (*this)[i] += scale * b[i];
                }
            }
        }
        virtual vec const& normalize() {
            return static_cast<vec&>(*this);
        }
        virtual base normalize(size_t i) {
            return (*this)[i];
        }
        void read() {
            for(auto &it: *this) {
                std::cin >> it;
            }
        }
        void print() const {
            std::ranges::copy(*this, std::ostream_iterator<base>(std::cout, " "));
            std::cout << "\n";
        }
        static vec random(size_t n) {
            vec res(n);
            std::ranges::generate(res, random::rng);
            return res;
        }
        // Concatenate vectors
        vec operator |(vec const& t) const {
            vec res(size(*this) + size(t));
            res[std::slice(0, size(*this), 1)] = *this;
            res[std::slice(size(*this), size(t), 1)] = t;
            return res;
        }

        // Generally, vec shouldn't be modified
        // after its pivot index is set
        std::pair<size_t, base> find_pivot() {
            if(pivot == size_t(-1)) {
                pivot = 0;
                while(pivot < size(*this) && normalize(pivot) == base(0)) {
                    pivot++;
                }
                if(pivot < size(*this)) {
                    pivot_inv = base(1) / (*this)[pivot];
                }
            }
            return {pivot, pivot_inv};
        }
        void reduce_by(vec &t) {
            auto [pivot, pinv] = t.find_pivot();
            if(pivot < size(*this)) {
                add_scaled(t, -normalize(pivot) * pinv, pivot);
            }
        }
    private:
        size_t pivot = -1;
        base pivot_inv;
    };

    template<typename base>
    struct vec: vec_base<vec<base>, base> {
        using Base = vec_base<vec<base>, base>;
        using Base::Base;
    };

    template<math::modint_type base>
    struct vec<base>: vec_base<vec<base>, base> {
        using Base = vec_base<vec<base>, base>;
        using Base::Base;

        void add_scaled(vec const& b, base scale, size_t i = 0) override {
            uint64_t scaler = scale.getr();
            if(scale != base(0)) {
                for(; i < size(*this); i++) {
                    (*this)[i].add_unsafe(scaler * b[i].getr_direct());
                }
                if(++counter == 8) {
                    for(auto &it: *this) {
                        it.pseudonormalize();
                    }
                    counter = 0;
                }
            }
        }
        vec const& normalize() override {
            for(auto &it: *this) {
                it.normalize();
            }
            return *this;
        }
        base normalize(size_t i) override {
            return (*this)[i].normalize();
        }
    private:
        size_t counter = 0;
    };
}
#endif // CP_ALGO_LINALG_VECTOR_HPP
