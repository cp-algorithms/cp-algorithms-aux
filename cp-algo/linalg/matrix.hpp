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
    template<typename base>
    struct matrix: std::valarray<vector<base>> {
        using Base = std::valarray<vector<base>>;
        using Base::Base;

        matrix(vector<base> t): Base(t, 1) {}
        matrix(size_t n): Base(vector<base>(n), n) {}
        matrix(size_t n, size_t m): Base(vector<base>(m), n) {}

        auto begin() {return std::begin(*static_cast<Base*>(this));}
        auto end() {return std::end(*static_cast<Base*>(this));}
        auto begin() const {return std::begin(*static_cast<Base const*>(this));}
        auto end() const {return std::end(*static_cast<Base const*>(this));}

        size_t n() const {return Base::size();}
        size_t m() const {return n() ? size(row(0)) : 0;}
        auto dim() const {return std::array{n(), m()};}

        auto& row(size_t i) {return (*this)[i];}
        auto const& row(size_t i) const {return (*this)[i];}

        matrix& operator *=(base t) {for(auto &it: *this) it *= t; return *this;}
        matrix operator *(base t) const {return matrix(*this) *= t;}

        matrix operator-() const {return Base::operator-();}
        matrix operator-(matrix const& t) const {return Base::operator-(t);}
        matrix operator+(matrix const& t) const {return Base::operator+(t);}
        matrix& operator*=(matrix const& t) {return *this = *this * t;}

        bool operator == (matrix const& t) const {
            return dim() == t.dim() && std::ranges::equal(*this, t,
                [](auto const& r1, auto const& r2) {
                    return (r1 == r2).min();
                });
        }
        bool operator != (matrix const& t) const {return !(*this == t);}

        void read() {
            std::ranges::for_each(*this, std::mem_fn(&vector<base>::read));
        }
        void print() const {
            std::ranges::for_each(*this, std::mem_fn(&vector<base>::print));
        }

        static matrix eye(size_t n) {
            matrix res(n);
            for(size_t i = 0; i < n; i++) {
                res[i][i] = 1;
            }
            return res;
        }

        // concatenate matrices
        matrix operator |(matrix const& b) const {
            assert(n() == b.n());
            matrix res(n(), m()+b.m());
            for(size_t i = 0; i < n(); i++) {
                res[i][std::slice(0, m(), 1)] = row(i);
                res[i][std::slice(m(), b.m(), 1)] = b[i];
            }
            return res;
        }
        matrix submatrix(auto slicex, auto slicey) const {
            matrix res = (*this)[slicex];
            for(auto &row: res) {
                row = vector<base>(row[slicey]);
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

        vector<base> apply(vector<base> const& x) const {
            return (matrix(x) * *this)[0];
        }

        matrix pow(uint64_t k) const {
            assert(n() == m());
            return bpow(*this, k, eye(n()));
        }

        static matrix rand(size_t n, size_t m) {
            matrix res(n, m);
            for(auto &it: res) {
                for(auto &jt: it) {
                    jt = random::rng();
                }
            }
            return res;
        }
        static matrix rand(size_t n) {
            return rand(n, n);
        }

        matrix& normalize() {
            for(auto &it: *this) {
                it.normalize();
            }
            return *this;
        }

        enum Mode {normal, reverse};
        template<Mode mode = normal>
        auto gauss(size_t lim) {
            size_t rk = 0;
            std::vector<size_t> free, pivots;
            for(size_t i = 0; i < lim; i++) {
                for(size_t j = rk; j < n() && row(rk).normalize(i) == base(0); j++) {
                    if(row(j).normalize(i) != 0) {
                        row(rk) += row(j);
                    }
                }
                if(rk == n() || row(rk).normalize()[i] == base(0)) {
                    free.push_back(i);
                } else {
                    pivots.push_back(i);
                    base dinv = -base(1) / row(rk)[i];
                    for(size_t j = (mode == normal) * rk; j < n(); j++) {
                        if(j != rk) {
                            row(j).add_scaled(row(rk), row(j).normalize(i) * dinv, i);
                        }
                    }
                    rk += 1;
                }
            }
            normalize();
            return std::array{pivots, free};
        }
        template<Mode mode = normal>
        auto gauss() {
            return gauss<mode>(m());
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
            b.gauss();
            base res = 1;
            for(size_t i = 0; i < n(); i++) {
                res *= b[i][i];
            }
            return res;
        }

        std::optional<matrix> inv() const {
            assert(n() == m());
            matrix b = *this | eye(n());
            if(size(b.gauss<reverse>(n())[0]) < n()) {
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
            auto [pivots, free] = A.template gauss<reverse>();
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
    };
}
#endif // CP_ALGO_LINALG_MATRIX_HPP
