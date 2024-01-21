
namespace algebra { // matrix
    template<int mod>
    struct matrix {
        using base = modular<mod>;
        size_t n, m;
        valarray<valarray<base>> a;
        matrix(size_t n, size_t m): n(n), m(m), a(valarray<base>(m), n) {}
        matrix(valarray<valarray<base>> a): n(size(a)), m(n ? size(a[0]) : 0), a(a) {}

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
                    cin >> (*this)[i][j];
                }
            }
        }

        void print() const {
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) {
                    cout << (*this)[i][j] << " \n"[j + 1 == m];
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
                res[i][slice(0,m,1)] = a[i];
                res[i][slice(m,b.m,1)] = b[i];
            }
            return res;
        }
        matrix submatrix(auto slicex, auto slicey) const {
            valarray res = a[slicex];
            for(auto &row: res) {
                row = valarray(row[slicey]);
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
                    for(size_t k = 0; k < b.m; k++) {
                        res[i][k].add_unsafe(a[i][j].r * b[j][k].r);
                    }
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

        inline static void mul_add(auto &a, auto const& b, base scale, size_t i = 0) {
            size_t m = size(a);
            for(; i < m; i++) {
                a[i].add_unsafe(scale.r * b[i].r);
            }
        }

        enum Mode {normal, reverse};
        template<Mode mode = normal>
        auto gauss(size_t lim) {
            size_t rk = 0;
            vector<size_t> free, pivots;
            for(size_t i = 0; i < lim; i++) {
                for(size_t j = rk; j < n && a[rk][i].normalize() == 0; j++) {
                    if(a[j][i].normalize() != 0) {
                        row(rk) += row(j);
                    }
                }
                if(rk == n || normalize(a[rk])[i] == 0) {
                    free.push_back(i);
                } else {
                    pivots.push_back(i);
                    base dinv = -a[rk][i].inv();
                    for(size_t j = mode == reverse ? 0 : rk; j < n; j++) {
                        if(j != rk) {
                            mul_add(a[j], a[rk], a[j][i].normalize() * dinv, i);
                        }
                    }
                    rk += 1;
                }
            }
            normalize();
            return array{pivots, free};
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

        optional<matrix> inv() const {
            assert(n == m);
            matrix b = *this | eye(n);
            if(size(b.gauss<reverse>(n)[0]) < n) {
                return nullopt;
            }
            for(size_t i = 0; i < n; i++) {
                b[i] *= b[i][i].inv();
            }
            return b.submatrix(slice(0, n, 1), slice(n, n, 1));
        }

        // [solution, basis], transposed
        optional<array<matrix, 2>> solve(matrix t) const {
            assert(n == t.n);
            matrix b = *this | t;
            auto [pivots, free] = b.gauss<reverse>();
            if(!empty(pivots) && pivots.back() >= m) {
                return nullopt;
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
            return array{
                sols.submatrix(slice(size(free) - t.m, t.m, 1), slice(0, m, 1)),
                sols.submatrix(slice(0, size(free) - t.m, 1), slice(0, m, 1))
            };
        }
    };
}
