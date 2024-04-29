#ifndef CP_ALGO_ALGEBRA_VECTOR_HPP
#define CP_ALGO_ALGEBRA_VECTOR_HPP
#include "common.hpp"
#include "modular.hpp"
#include <functional>
#include <algorithm>
#include <valarray>
#include <iostream>
#include <iterator>
namespace cp_algo::algebra {
    template<class derived, typename base>
    struct vector_base: std::valarray<base> {
        using Base = std::valarray<base>;
        using Base::Base;

        auto begin() {return std::begin(*static_cast<Base*>(this));}
        auto end() {return std::end(*static_cast<Base*>(this));}
        auto begin() const {return std::begin(*static_cast<Base const*>(this));}
        auto end() const {return std::end(*static_cast<Base const*>(this));}

        void add_scaled(derived const& b, base scale, size_t i = 0) {
            for(; i < this->size(); i++) {
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
    };

    template<typename base>
    struct vector: vector_base<vector<base>, base> {
        using Base = vector_base<vector<base>, base>;
        using Base::Base;
    };

    template<int mod>
    struct vector<modular<mod>>: vector_base<vector<modular<mod>>, modular<mod>> {
        using base = modular<mod>;
        using Base = vector_base<vector<base>, base>;
        using Base::Base;

        size_t counter = 0;
        void add_scaled(vector const& b, base scale, size_t i = 0) {
            for(; i < this->size(); i++) {
                (*this)[i].add_unsafe(scale.r * b[i].r);
            }
            counter++;
            if(counter == 8) {
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
    };
}
#endif // CP_ALGO_ALGEBRA_VECTOR_HPP
