#ifndef CP_ALGO_LINALG_VECTOR_HPP
#define CP_ALGO_LINALG_VECTOR_HPP
#include "../algebra/modular.hpp"
#include "../random/rng.hpp"
#include <functional>
#include <algorithm>
#include <valarray>
#include <iostream>
#include <iterator>
namespace cp_algo::linalg {
    template<class derived, typename base>
    struct valarray_base: std::valarray<base> {
        using Base = std::valarray<base>;
        using Base::Base;

        valarray_base(base const& t): Base(t, 1) {}

        auto begin() {return std::begin(*static_cast<Base*>(this));}
        auto end() {return std::end(*static_cast<Base*>(this));}
        auto begin() const {return std::begin(*static_cast<Base const*>(this));}
        auto end() const {return std::end(*static_cast<Base const*>(this));}

        bool operator == (derived const& t) const {return std::ranges::equal(*this, t);}
        bool operator != (derived const& t) const {return !(*this == t);}

        derived operator-() const {return Base::operator-();}
        derived operator-(derived const& t) const {return Base::operator-(t);}
        derived operator+(derived const& t) const {return Base::operator+(t);}

        static derived from_range(auto const& R) {
            derived res(std::ranges::distance(R));
            std::ranges::copy(R, res.begin());
            return res;
        }
    };

    template<class vector, typename base>
    struct vector_base: valarray_base<vector_base<vector, base>, base> {
        using Base = valarray_base<vector_base<vector, base>, base>;
        using Base::Base;

        static vector ei(size_t n, size_t i) {
            vector res(n);
            res[i] = 1;
            return res;
        }

        // Make sure the result is vector, not Base
        vector operator*(base t) const {return Base::operator*(t);}

        void add_scaled(vector const& b, base scale, size_t i = 0) {
            for(; i < size(*this); i++) {
                (*this)[i] += scale * b[i];
            }
        }
        auto& normalize() {
            return *this;
        }
        auto& normalize(size_t i) {
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
        static vector random(size_t n) {
            vector res(n);
            std::ranges::generate(res, random::rng);
            return res;
        }
        // Concatenate vectors
        vector operator |(vector const& t) const {
            vector res(size(*this) + size(t));
            res[std::slice(0, size(*this), 1)] = *this;
            res[std::slice(size(*this), size(t), 1)] = t;
            return res;
        }

        // Generally, vector shouldn't be modified
        // after it's pivot index is set
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
        void reduce_by(vector &t) {
            auto [pivot, pinv] = t.find_pivot();
            if(pivot < size(*this)) {
                add_scaled(t, -(*this)[pivot] * pinv, pivot);
            }
        }
    private:
        size_t pivot = -1;
        base pivot_inv;
    };

    template<typename base>
    struct vec: vector_base<vec<base>, base> {
        using Base = vector_base<vec<base>, base>;
        using Base::Base;
    };

    template<int mod>
    struct vec<algebra::modular<mod>>:
            vector_base<vec<algebra::modular<mod>>, algebra::modular<mod>> {
        using base = algebra::modular<mod>;
        using Base = vector_base<vec<base>, base>;
        using Base::Base;

        void add_scaled(vec const& b, base scale, size_t i = 0) {
            for(; i < size(*this); i++) {
                (*this)[i].add_unsafe(scale.r * b[i].r);
            }
            if(++counter == 8) {
                std::ranges::for_each(*this, std::mem_fn(&base::pseudonormalize));
                counter = 0;
            }
        }
        auto& normalize() {
            std::ranges::for_each(*this, std::mem_fn(&base::normalize));
            return *this;
        }
        auto& normalize(size_t i) {
            return (*this)[i].normalize();
        }
    private:
        size_t counter = 0;
    };
}
#endif // CP_ALGO_LINALG_VECTOR_HPP
