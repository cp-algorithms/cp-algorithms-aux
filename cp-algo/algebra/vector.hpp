#ifndef CP_ALGO_ALGEBRA_VECTOR_HPP
#define CP_ALGO_ALGEBRA_VECTOR_HPP
#include "common.hpp"
#include "modular.hpp"
#include <valarray>
namespace cp_algo::algebra {
    template<class derived, typename base>
    struct vector_base: std::valarray<base> {
        using Base = std::valarray<base>;
        using Base::Base;

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
                counter = 0;
                for(size_t i = 0; i < this->size(); i++) {
                    (*this)[i].pseudonormalize();
                }
            }
        }
        auto& normalize() {
            for(size_t i = 0; i < this->size(); i++) {
                (*this)[i].normalize();
            }
            return *this;
        }
        auto& normalize(size_t i) {
            return (*this)[i].normalize();
        }
    };
}
#endif // CP_ALGO_ALGEBRA_VECTOR_HPP
