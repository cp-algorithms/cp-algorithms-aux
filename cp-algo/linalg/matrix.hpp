#ifndef CP_ALGO_LINALG_MATRIX_HPP
#define CP_ALGO_LINALG_MATRIX_HPP
#include "../random/rng.hpp"
#include "../algebra/common.hpp"
#include "vector.hpp"
#include <iostream>
#include <optional>
#include <cassert>
#include <vector>
#include <array>
namespace cp_algo::linalg {
    template<typename base_t>
    struct matrix: valarray_base<matrix<base_t>, vec<base_t>> {
        using base = base_t;
        using Base = valarray_base<matrix<base>, vec<base>>;
        using Base::Base;

        matrix(size_t n): Base(vec<base>(n), n) {}
        matrix(size_t n, size_t m): Base(vec<base>(m), n) {}

        size_t n() const {return size(*this);}
        size_t m() const {return n() ? size(row(0)) : 0;}
        auto dim() const {return std::array{n(), m()};}

        auto& row(size_t i) {return (*this)[i];}
        auto const& row(size_t i) const {return (*this)[i];}

        matrix& operator *=(base t) {for(auto &it: *this) it *= t; return *this;}
        matrix operator *(base t) const {return matrix(*this) *= t;}

        // Make sure the result is matrix, not Base
        matrix& operator*=(matrix const& t) {return *this = *this * t;}

        void read() {
            for(auto &it: *this) {
                it.read();
            }
        }
        void print() const {
            for(auto const& it: *this) {
                it.print();
            }
        }

        static matrix block_diagonal(std::vector<matrix> const& blocks) {
            size_t n = 0;
            for(auto &it: blocks) {
                assert(it.n() == it.m());
                n += it.n();
            }
            matrix res(n);
            n = 0;
            for(auto &it: blocks) {
                for(size_t i = 0; i < it.n(); i++) {
                    res[n + i][std::slice(n, it.n(), 1)] = it[i];
                }
                n += it.n();
            }
            return res;
        }
        static matrix random(size_t n, size_t m) {
            matrix res(n, m);
            std::ranges::generate(res, std::bind(vec<base>::random, m));
            return res;
        }
        static matrix random(size_t n) {
            return random(n, n);
        }
        static matrix eye(size_t n) {
            matrix res(n);
            for(size_t i = 0; i < n; i++) {
                res[i][i] = 1;
            }
            return res;
        }

        // Concatenate matrices
        matrix operator |(matrix const& b) const {
            assert(n() == b.n());
            matrix res(n(), m()+b.m());
            for(size_t i = 0; i < n(); i++) {
                res[i] = row(i) | b[i];
            }
            return res;
        }
        matrix submatrix(auto slicex, auto slicey) const {
            matrix res = (*this)[slicex];
            for(auto &row: res) {
                row = vec<base>(row[slicey]);
            }
            return res;
        }

        matrix T() const {
            matrix res(m(), n());
            for(size_t i = 0; i < n(); i++) {
                for(size_t j = 0; j < m(); j++) {
                    res[j][i] = row(i)[j];
                }
            }
            return res;
        }

        matrix operator *(matrix const& b) const {
            assert(m() == b.n());
            matrix res(n(), b.m());
            for(size_t i = 0; i < n(); i++) {
                for(size_t j = 0; j < m(); j++) {
                    res[i].add_scaled(b[j], row(i)[j]);
                }
            }
            return res.normalize();
        }

        vec<base> apply(vec<base> const& x) const {
            return (matrix(x) * *this)[0];
        }

        matrix pow(uint64_t k) const {
            assert(n() == m());
            return bpow(*this, k, eye(n()));
        }

        matrix& normalize() {
            for(auto &it: *this) {
                it.normalize();
            }
            return *this;
        }

        enum gauss_mode {normal, reverse};
        template<gauss_mode mode = normal>
        matrix& gauss() {
            for(size_t i = 0; i < n(); i++) {
                row(i).normalize();
                for(size_t j = (mode == normal) * i; j < n(); j++) {
                    if(j != i) {
                        row(j).reduce_by(row(i));
                    }
                }
            }
            return normalize();
        }
        template<gauss_mode mode = normal>
        auto echelonize(size_t lim) {
            return gauss<mode>().sort_classify(lim);
        }
        template<gauss_mode mode = normal>
        auto echelonize() {
            return echelonize<mode>(m());
        }

        size_t rank() const {
            if(n() < m()) {
                return T().rank();
            }
            return size(matrix(*this).gauss()[0]);
        }

        base det() const {
            assert(n() == m());
            matrix b = *this;
            b.echelonize();
            base res = 1;
            for(size_t i = 0; i < n(); i++) {
                res *= b[i][i];
            }
            return res;
        }

        std::optional<matrix> inv() const {
            assert(n() == m());
            matrix b = *this | eye(n());
            if(size(b.echelonize<reverse>(n())[0]) < n()) {
                return std::nullopt;
            }
            for(size_t i = 0; i < n(); i++) {
                b[i] *= base(1) / b[i][i];
            }
            return b.submatrix(std::slice(0, n(), 1), std::slice(n(), n(), 1));
        }

        // Can also just run gauss on T() | eye(m)
        // but it would be slower :(
        auto kernel() const {
            auto A = *this;
            auto [pivots, free] = A.template echelonize<reverse>();
            matrix sols(size(free), m());
            for(size_t j = 0; j < size(pivots); j++) {
                base scale = base(1) / A[j][pivots[j]];
                for(size_t i = 0; i < size(free); i++) {
                    sols[i][pivots[j]] = A[j][free[i]] * scale;
                }
            }
            for(size_t i = 0; i < size(free); i++) {
                sols[i][free[i]] = -1;
            }
            return sols;
        }

        // [solution, basis], transposed
        std::optional<std::array<matrix, 2>> solve(matrix t) const {
            matrix sols = (*this | t).kernel();
            if(sols.n() < t.m() || sols.submatrix(
                std::slice(sols.n() - t.m(), t.m(), 1),
                std::slice(m(), t.m(), 1)
            ) != -eye(t.m())) {
                return std::nullopt;
            } else {
                return std::array{
                    sols.submatrix(std::slice(sols.n() - t.m(), t.m(), 1),
                                   std::slice(0, m(), 1)),
                    sols.submatrix(std::slice(0, sols.n() - t.m(), 1),
                                   std::slice(0, m(), 1))
                };
            }
        }
    private:
        // To be called after a gaussian elimination run
        // Sorts rows by pivots and classifies
        // variables into pivots and free
        auto sort_classify(size_t lim) {
            size_t rk = 0;
            std::vector<size_t> free, pivots;
            for(size_t j = 0; j < lim; j++) {
                for(size_t i = rk + 1; i < n() && row(rk)[j] == base(0); i++) {
                    if(row(i)[j] != base(0)) {
                        std::swap(row(i), row(rk));
                        row(rk) = -row(rk);
                    }
                }
                if(rk < n() && row(rk)[j] != base(0)) {
                    pivots.push_back(j);
                    rk++;
                } else {
                    free.push_back(j);
                }
            }
            return std::array{pivots, free};
        }
    };
}
#endif // CP_ALGO_LINALG_MATRIX_HPP
