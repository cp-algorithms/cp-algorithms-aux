#ifndef CP_ALGO_ALGEBRA_MATRIX_HPP
#define CP_ALGO_ALGEBRA_MATRIX_HPP
#include "common.hpp"
#include "modular.hpp"
#include <valarray>
#include <iostream>
#include <optional>
#include <cassert>
#include <vector>
#include <array>
namespace cp_algo::algebra {
    template<int mod>
    struct matrix {
        using base = modular<mod>;
        size_t n, m;
        std::valarray<std::valarray<base>> a;
        matrix(size_t n, size_t m): n(n), m(m), a(std::valarray<base>(m), n) {}
        matrix(std::valarray<std::valarray<base>> a): n(size(a)), m(n ? size(a[0]) : 0), a(a) {}

        auto& operator[] (size_t i) {return a[i];}
        auto const& operator[] (size_t i) const {return a[i];}
        auto& row(size_t i) {return a[i];}
        auto const& row(size_t i) const {return a[i];}

        matrix operator -() const {return matrix(-a);}
        matrix& operator *=(base t) {for(auto &it: a) it *= t; return *this;}
        matrix operator *(base t) const {return matrix(*this) *= t;}

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
                row = std::valarray(row[slicey]);
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
                    add_scaled(res[i], b[j], a[i][j]);
                }
            }
            return res.normalize();
        }

        matrix pow(uint64_t k) const {
            assert(n == m);
            return bpow(*this, k, eye(n));
        }

        static auto& normalize(auto &a) {
            for(auto &it: a) {
                it.normalize();
            }
            return a;
        }
        matrix& normalize() {
            for(auto &it: a) {
                normalize(it);
            }
            return *this;
        }

        inline static void add_scaled(auto &a, auto const& b, base scale, size_t i = 0) {
            size_t m = size(a);
            for(; i < m; i++) {
                a[i].add_unsafe(scale.r * b[i].r);
            }
        }

        enum Mode {normal, reverse};
        template<Mode mode = normal>
        auto gauss(size_t lim) {
            size_t rk = 0;
            std::vector<size_t> free, pivots;
            for(size_t i = 0; i < lim; i++) {
                for(size_t j = rk; j < n && a[rk][i].normalize() == 0; j++) {
                    if(a[j][i].normalize() != 0) {
                        a[rk] += a[j];
                    }
                }
                if(rk == n || normalize(a[rk])[i] == 0) {
                    free.push_back(i);
                } else {
                    pivots.push_back(i);
                    base dinv = -a[rk][i].inv();
                    for(size_t j = mode == reverse ? 0 : rk; j < n; j++) {
                        if(j != rk) {
                            add_scaled(a[j], a[rk], a[j][i].normalize() * dinv, i);
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
                b[i] *= b[i][i].inv();
            }
            return b.submatrix(std::slice(0, n, 1), std::slice(n, n, 1));
        }

        // [solution, basis], transposed
        std::optional<std::array<matrix, 2>> solve(matrix t) const {
            assert(n == t.n);
            matrix b = *this | t;
            auto [pivots, free] = b.gauss<reverse>();
            if(!empty(pivots) && pivots.back() >= m) {
                return std::nullopt;
            }
            matrix sols(size(free), m);
            for(size_t j = 0; j < size(pivots); j++) {
                base scale = b[j][pivots[j]].inv();
                for(size_t i = 0; i < size(free); i++) {
                    sols[i][pivots[j]] = b[j][free[i]] * scale;
                }
            }
            for(size_t i = 0; free[i] < m; i++) {
                sols[i][free[i]] = -1;
            }
            return std::array{
                sols.submatrix(std::slice(size(free) - t.m, t.m, 1), std::slice(0, m, 1)),
                sols.submatrix(std::slice(0, size(free) - t.m, 1), std::slice(0, m, 1))
            };
        }
    };
}
#endif // CP_ALGO_ALGEBRA_MATRIX_HPP
