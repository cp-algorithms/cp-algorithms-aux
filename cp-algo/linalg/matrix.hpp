#ifndef CP_ALGO_LINALG_MATRIX_HPP
#define CP_ALGO_LINALG_MATRIX_HPP
#include "../random/rng.hpp"
#include "../math/common.hpp"
#include "vector.hpp"
#include <iostream>
#include <optional>
#include <cassert>
#include <vector>
#include <array>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::linalg {
    enum gauss_mode {normal, reverse};

    template<typename base_t, class _vec_t = std::conditional_t<
        math::modint_type<base_t>,
        modint_vec<base_t>,
        vec<base_t>>>
    struct matrix: big_vector<_vec_t> {
        using vec_t = _vec_t;
        using base = base_t;
        using Base = big_vector<vec_t>;
        using Base::Base;

        matrix(size_t n): Base(n, vec_t(n)) {}
        matrix(size_t n, size_t m): Base(n, vec_t(m)) {}

        matrix(Base const& t): Base(t) {}
        matrix(Base &&t): Base(std::move(t)) {}
        
        template<std::ranges::input_range R>
        matrix(R &&r): Base(r.begin(), r.end()) {}

        size_t n() const {return size(*this);}
        size_t m() const {return n() ? size(row(0)) : 0;}
        
        void resize(size_t n, size_t m) {
            Base::resize(n);
            for(auto &it: *this) {
                it.resize(m);
            }
        }

        auto& row(size_t i) {return (*this)[i];}
        auto const& row(size_t i) const {return (*this)[i];}

        auto elements() {return *this | std::views::join;}
        auto elements() const {return *this | std::views::join;}

        matrix operator-() const {
            return *this | std::views::transform([](auto x) {return vec_t(-x);});
        }
        matrix& operator+=(matrix const& t) {
            for(auto [a, b]: std::views::zip(elements(), t.elements())) {
                a += b;
            }
            return *this;
        }
        matrix& operator -=(matrix const& t) {
            for(auto [a, b]: std::views::zip(elements(), t.elements())) {
                a -= b;
            }
            return *this;
        }
        matrix operator+(matrix const& t) const {return matrix(*this) += t;}
        matrix operator-(matrix const& t) const {return matrix(*this) -= t;}
        
        matrix& operator *=(base t) {for(auto &it: *this) it *= t; return *this;}
        matrix operator *(base t) const {return matrix(*this) *= t;}
        matrix& operator /=(base t) {return *this *= base(1) / t;}
        matrix operator /(base t) const {return matrix(*this) /= t;}

        // Make sure the result is matrix, not Base
        matrix& operator *=(matrix const& t) {return *this = *this * t;}

        void read_transposed() {
            for(size_t j = 0; j < m(); j++) {
                for(size_t i = 0; i < n(); i++) {
                    std::cin >> (*this)[i][j];
                }
            }
        }
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

        static matrix block_diagonal(big_vector<matrix> const& blocks) {
            size_t n = 0;
            for(auto &it: blocks) {
                assert(it.n() == it.m());
                n += it.n();
            }
            matrix res(n);
            n = 0;
            for(auto &it: blocks) {
                for(size_t i = 0; i < it.n(); i++) {
                    std::ranges::copy(it[i], begin(res[n + i]) + n);
                }
                n += it.n();
            }
            return res;
        }
        static matrix random(size_t n, size_t m) {
            matrix res(n, m);
            std::ranges::generate(res, std::bind(vec_t::random, m));
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
        void assign_submatrix(auto viewx, auto viewy, matrix const& t) {
            for(auto [a, b]: std::views::zip(*this | viewx, t)) {
                std::ranges::copy(b, begin(a | viewy));
            }
        }
        auto submatrix(auto viewx, auto viewy) const {
            return *this | viewx | std::views::transform([viewy](auto const& y) {
                return y | viewy;
            });
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

        vec_t apply(vec_t const& x) const {
            return (matrix(1, x) * *this)[0];
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
        template<gauss_mode mode = normal>
        void eliminate(size_t i, size_t k) {
            auto kinv = base(1) / row(i).normalize()[k];
            for(size_t j = (mode == normal) * i; j < n(); j++) {
                if(j != i) {
                    row(j).add_scaled(row(i), -row(j).normalize(k) * kinv);
                }
            }
        }
        template<gauss_mode mode = normal>
        void eliminate(size_t i) {
            row(i).normalize();
            for(size_t j = (mode == normal) * i; j < n(); j++) {
                if(j != i) {
                    row(j).reduce_by(row(i));
                }
            }
        }
        template<gauss_mode mode = normal>
        matrix& gauss() {
            for(size_t i = 0; i < n(); i++) {
                eliminate<mode>(i);
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
            if(n() > m()) {
                return T().rank();
            }
            return size(matrix(*this).echelonize()[0]);
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

        std::pair<base, matrix> inv() const {
            assert(n() == m());
            matrix b = *this | eye(n());
            if(size(b.echelonize<reverse>(n())[0]) < n()) {
                return {0, {}};
            }
            base det = 1;
            for(size_t i = 0; i < n(); i++) {
                det *= b[i][i];
                b[i] *= base(1) / b[i][i];
            }
            return {det, b.submatrix(std::views::all, std::views::drop(n()))};
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
            if(sols.n() < t.m() || matrix(sols.submatrix(
                std::views::drop(sols.n() - t.m()),
                std::views::drop(m())
            )) != -eye(t.m())) {
                return std::nullopt;
            } else {
                return std::array{
                    matrix(sols.submatrix(std::views::drop(sols.n() - t.m()), std::views::take(m()))),
                    matrix(sols.submatrix(std::views::take(sols.n() - t.m()), std::views::take(m())))
                };
            }
        }

        // To be called after a gaussian elimination run
        // Sorts rows by pivots and classifies
        // variables into pivots and free
        auto sort_classify(size_t lim) {
            size_t rk = 0;
            big_vector<size_t> free, pivots;
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
    template<typename base_t>
    auto operator *(base_t t, matrix<base_t> const& A) {return A * t;}
}
#pragma GCC pop_options
#endif // CP_ALGO_LINALG_MATRIX_HPP
