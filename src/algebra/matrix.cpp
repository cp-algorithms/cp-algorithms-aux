namespace algebra { // matrix
    template<int mod>
    struct matrix {
        using _ = modular<mod>;
        size_t n, m;
        valarray<valarray<_>> a;
        matrix(size_t n, size_t m): n(n), m(m), a(valarray<_>(m), n) {}
        matrix(valarray<valarray<_>> a): n(size(a)), m(n ? size(a[0]) : 0), a(a) {}

        auto& operator[] (size_t i) {return a[i];}
        auto const& operator[] (size_t i) const {return a[i];}
        auto& row(size_t i) {return a[i];}
        auto const& row(size_t i) const {return a[i];}

        matrix operator -() const {return matrix(-a);}
        matrix& operator *=(_ t) {for(auto &it: a) it *= t; return *this;}
        matrix operator *(_ t) const {return matrix(*this) *= t;}

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
        matrix submatrix(auto slicex, auto slicey) {
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
                        res[i][k].r += a[i][j].r * b[j][k].r;
                        if((m - j) % 16 == 1) {
                            res[i][k].r %= mod;
                        }
                    }
                }
            }
            return res;
        }

        matrix pow(uint64_t k) const {
            assert(n == m);
            return bpow(*this, k, eye(n));
        }

        void normalize(size_t i) {
            for(size_t j = 0; j < m; j++) {
                a[i][j].r %= mod;
            }
        }

        template<bool reverse = false>
        auto gauss(size_t lim) {
            vector<size_t> free, pivots;
            size_t rk = 0;
            for(size_t i = 0; i < lim; i++) {
                for(size_t j = rk; j < n; j++) {
                    a[j][i].r %= mod;
                    if(a[j][i] != 0) {
                        normalize(j);
                        if(rk != j) {
                            swap(a[rk], a[j]);
                            a[rk] *= -1;
                        }
                        break;
                    }
                }
                if(rk == n || a[rk][i] == 0) {
                    free.push_back(i);
                    continue;
                }
                _ dinv = -a[rk][i].inv();
                for(size_t j = reverse ? 0 : rk + 1; j < n; j++) {
                    if(j != rk) {
                        a[j][i].r %= mod;
                        _ scale = a[j][i] * dinv;
                        for(size_t k = i; k < m; k++) {
                            a[j][k].r += scale.r * a[rk][k].r;
                            if(rk % 16 == 15) {
                                a[j][k].r %= mod;
                            }
                        }
                    }
                }
                pivots.push_back(i);
                rk += 1;
            }
            for(size_t i = 0; i < n; i++) {
                normalize(i);
            }
            return array{pivots, free};
        }
        template<bool reverse = false>
        auto gauss() {
            return gauss<reverse>(m);
        }

        size_t rank() const {
            if(n < m) {
                return T().rank();
            }
            return size(matrix(*this).gauss()[0]);
        }

        _ det() const {
            assert(n == m);
            matrix b = *this;
            b.gauss();
            _ res = 1;
            for(size_t i = 0; i < n; i++) {
                res *= b[i][i];
            }
            return res;
        }

        optional<matrix> inv() const {
            assert(n == m);
            matrix b = *this | eye(n);
            if(size(b.gauss<true>(n)[0]) < n) {
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
            auto [pivots, free] = b.gauss<true>();
            if(!empty(pivots) && pivots.back() >= m) {
                return nullopt;
            }
            matrix sols(size(free), m);
            for(size_t j = 0; j < size(pivots); j++) {
                _ scale = b[j][pivots[j]].inv();
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
