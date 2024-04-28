#ifndef CP_ALGO_ALGEBRA_MATRIX_HPP
#define CP_ALGO_ALGEBRA_MATRIX_HPP
#include "../random/rng.hpp"
#include "vector.hpp"
#include <iostream>
#include <optional>
#include <cassert>
#include <vector>
#include <array>
namespace cp_algo::algebra {
    template<typename base>
    struct matrix {
        size_t n, m;
        std::valarray<vector<base>> a;
        matrix(vector<base> t): n(1), m(size(t)), a(t, 1) {}
        matrix(size_t n, size_t m): n(n), m(m), a(vector<base>(m), n) {}
        matrix(std::valarray<vector<base>> a): n(size(a)), m(n ? size(a[0]) : 0), a(a) {}

        auto& operator[] (size_t i) {return a[i];}
        auto const& operator[] (size_t i) const {return a[i];}
        auto& row(size_t i) {return a[i];}
        auto const& row(size_t i) const {return a[i];}

        matrix operator -() const {return matrix(-a);}
        matrix& operator *=(base t) {for(auto &it: a) it *= t; return *this;}
        matrix operator *(base t) const {return matrix(*this) *= t;}

        bool operator == (matrix const& t) const {
            if(n != t.n || m != t.m) {
                return false;
            }
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) {
                    if(a[i][j] != t.a[i][j]) {
                        return false;
                    }
                }
            }
            return true;
        }
        bool operator != (matrix const& t) const {return !(*this == t);}

        void read() {
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) {
                    std::cin >> (*this)[i][j];
                }
            }
        }

        void print() const {
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) {
                    std::cout << (*this)[i][j] << " \n"[j + 1 == m];
                }
            }
        }

        static matrix eye(size_t n) {
            matrix res(n, n);
            for(size_t i = 0; i < n; i++) {
                res[i][i] = 1;
            }
            return res;
        }

        // concatenate matrices
        matrix operator |(matrix const& b) const {
            assert(n == b.n);
            matrix res(n, m+b.m);
            for(size_t i = 0; i < n; i++) {
                res[i][std::slice(0,m,1)] = a[i];
                res[i][std::slice(m,b.m,1)] = b[i];
            }
            return res;
        }
        matrix submatrix(auto slicex, auto slicey) const {
            std::valarray res = a[slicex];
            for(auto &row: res) {
                row = vector<base>(row[slicey]);
            }
            return res;
        }

        matrix T() const {
            matrix res(m, n);
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) {
                    res[j][i] = (*this)[i][j];
                }
            }
            return res;
        }

        matrix operator *(matrix const& b) const {
            assert(m == b.n);
            matrix res(n, b.m);
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) {
                    res[i].add_scaled(b[j], a[i][j]);
                }
            }
            return res.normalize();
        }

        vector<base> apply(vector<base> const& x) const {
            return (matrix(x) * *this)[0];
        }

        matrix pow(uint64_t k) const {
            assert(n == m);
            return bpow(*this, k, eye(n));
        }

        static matrix rand(size_t n, size_t m) {
            matrix res(n, m);
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) {
                    res[i][j] = random::rng();
                }
            }
            return res;
        }
        static vector<base> rand(size_t n) {
            return rand(1, n)[0];
        }

        matrix& normalize() {
            for(auto &it: a) {
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
                for(size_t j = rk; j < n && a[rk].normalize(i) == 0; j++) {
                    if(a[j].normalize(i) != 0) {
                        a[rk] += a[j];
                    }
                }
                if(rk == n || a[rk].normalize()[i] == 0) {
                    free.push_back(i);
                } else {
                    pivots.push_back(i);
                    base dinv = -base(1) / a[rk][i];
                    for(size_t j = mode == reverse ? 0 : rk; j < n; j++) {
                        if(j != rk) {
                            a[j].add_scaled(a[rk], a[j].normalize(i) * dinv, i);
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
            return gauss<mode>(m);
        }

        size_t rank() const {
            if(n < m) {
                return T().rank();
            }
            return size(matrix(*this).gauss()[0]);
        }

        base det() const {
            assert(n == m);
            matrix b = *this;
            b.gauss();
            base res = 1;
            for(size_t i = 0; i < n; i++) {
                res *= b[i][i];
            }
            return res;
        }

        std::optional<matrix> inv() const {
            assert(n == m);
            matrix b = *this | eye(n);
            if(size(b.gauss<reverse>(n)[0]) < n) {
                return std::nullopt;
            }
            for(size_t i = 0; i < n; i++) {
                b[i] *= base(1) / b[i][i];
            }
            return b.submatrix(std::slice(0, n, 1), std::slice(n, n, 1));
        }

        // Can also just run gauss on T() | eye(m)
        // but it would be slower :(
        auto kernel() const {
            auto A = *this;
            auto [pivots, free] = A.template gauss<reverse>();
            matrix sols(size(free), m);
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
            if(sols.n < t.m || sols.submatrix(
                std::slice(sols.n - t.m, t.m, 1),
                std::slice(m, t.m, 1)
            ) != -eye(t.m)) {
                return std::nullopt;
            } else {
                return std::array{
                    sols.submatrix(std::slice(sols.n - t.m, t.m, 1),
                                   std::slice(0, m, 1)),
                    sols.submatrix(std::slice(0, sols.n - t.m, 1),
                                   std::slice(0, m, 1))
                };
            }
        }
    };
}
#endif // CP_ALGO_ALGEBRA_MATRIX_HPP
