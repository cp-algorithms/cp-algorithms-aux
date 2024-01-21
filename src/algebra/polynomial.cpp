namespace algebra { // poly
    template<typename T>
    struct poly {
        vector<T> a;
        
        void normalize() { // get rid of leading zeroes
            while(!a.empty() && a.back() == T(0)) {
                a.pop_back();
            }
        }
        
        poly(){}
        poly(T a0) : a{a0}{normalize();}
        poly(const vector<T> &t) : a(t){normalize();}
        
        poly operator -() const {
            auto t = *this;
            for(auto &it: t.a) {
                it = -it;
            }
            return t;
        }
        
        poly operator += (const poly &t) {
            a.resize(max(a.size(), t.a.size()));
            for(size_t i = 0; i < t.a.size(); i++) {
                a[i] += t.a[i];
            }
            normalize();
            return *this;
        }
        
        poly operator -= (const poly &t) {
            a.resize(max(a.size(), t.a.size()));
            for(size_t i = 0; i < t.a.size(); i++) {
                a[i] -= t.a[i];
            }
            normalize();
            return *this;
        }
        poly operator + (const poly &t) const {return poly(*this) += t;}
        poly operator - (const poly &t) const {return poly(*this) -= t;}
        
        poly mod_xk(size_t k) const { // get first k coefficients
            return vector<T>(begin(a), begin(a) + min(k, a.size()));
        }
        
        poly mul_xk(size_t k) const { // multiply by x^k
            auto res = a;
            res.insert(begin(res), k, 0);
            return res;
        }
        
        poly div_xk(size_t k) const { // drop first k coefficients
            return vector<T>(begin(a) + min(k, a.size()), end(a));
        }
        
        poly substr(size_t l, size_t r) const { // return mod_xk(r).div_xk(l)
            return vector<T>(
                begin(a) + min(l, a.size()),
                begin(a) + min(r, a.size())
            );
        }
        
        poly operator *= (const poly &t) {fft::mul(a, t.a); normalize(); return *this;}
        poly operator * (const poly &t) const {return poly(*this) *= t;}
        
        poly reverse(size_t n) const { // computes x^n A(x^{-1})
            auto res = a;
            res.resize(max(n, res.size()));
            return vector<T>(res.rbegin(), res.rbegin() + n);
        }
        
        poly reverse() const {
            return reverse(deg() + 1);
        }
        
        pair<poly, poly> divmod_slow(const poly &b) const { // when divisor or quotient is small
            vector<T> A(a);
            vector<T> res;
            T b_lead_inv = b.a.back().inv();
            while(A.size() >= b.a.size()) {
                res.push_back(A.back() * b_lead_inv);
                if(res.back() != T(0)) {
                    for(size_t i = 0; i < b.a.size(); i++) {
                        A[A.size() - i - 1] -= res.back() * b.a[b.a.size() - i - 1];
                    }
                }
                A.pop_back();
            }
            std::reverse(begin(res), end(res));
            return {res, A};
        }
        
        pair<poly, poly> divmod_hint(poly const& b, poly const& binv) const { // when inverse is known
            assert(!b.is_zero());
            if(deg() < b.deg()) {
                return {poly{0}, *this};
            }
            int d = deg() - b.deg();
            if(min(d, b.deg()) < magic) {
                return divmod_slow(b);
            }
            poly D = (reverse().mod_xk(d + 1) * binv.mod_xk(d + 1)).mod_xk(d + 1).reverse(d + 1);
            return {D, *this - D * b};
        }
        
        pair<poly, poly> divmod(const poly &b) const { // returns quotiend and remainder of a mod b
            assert(!b.is_zero());
            if(deg() < b.deg()) {
                return {poly{0}, *this};
            }
            int d = deg() - b.deg();
            if(min(d, b.deg()) < magic) {
                return divmod_slow(b);
            }
            poly D = (reverse().mod_xk(d + 1) * b.reverse().inv(d + 1)).mod_xk(d + 1).reverse(d + 1);
            return {D, *this - D * b};
        }
        
        // (ax+b) / (cx+d)
        struct transform {
            poly a, b, c, d;
            transform(poly a, poly b = T(1), poly c = T(1), poly d = T(0)): a(a), b(b), c(c), d(d){}
            
            transform operator *(transform const& t) {
                return {
                    a*t.a + b*t.c, a*t.b + b*t.d,
                    c*t.a + d*t.c, c*t.b + d*t.d
                };
            }
            
            transform adj() {
                return transform(d, -b, -c, a);
            }
            
            auto apply(poly A, poly B) {
                return make_pair(a * A + b * B, c * A + d * B);
            }
        };
        
        template<typename Q>
        static void concat(vector<Q> &a, vector<Q> const& b) {
            for(auto it: b) {
                a.push_back(it);
            }
        }
        
        // finds a transform that changes A/B to A'/B' such that
        // deg B' is at least 2 times less than deg A
        static pair<vector<poly>, transform> half_gcd(poly A, poly B) {
            assert(A.deg() >= B.deg());
            int m = (A.deg() + 1) / 2;
            if(B.deg() < m) {
                return {{}, {T(1), T(0), T(0), T(1)}};
            }
            auto [ar, Tr] = half_gcd(A.div_xk(m), B.div_xk(m));
            tie(A, B) = Tr.adj().apply(A, B);
            if(B.deg() < m) {
                return {ar, Tr};
            }
            auto [ai, R] = A.divmod(B);
            tie(A, B) = make_pair(B, R);
            int k = 2 * m - B.deg();
            auto [as, Ts] = half_gcd(A.div_xk(k), B.div_xk(k));
            concat(ar, {ai});
            concat(ar, as);
            return {ar, Tr * transform(ai) * Ts};
        }
        
        // return a transform that reduces A / B to gcd(A, B) / 0
        static pair<vector<poly>, transform> full_gcd(poly A, poly B) {
            vector<poly> ak;
            vector<transform> trs;
            while(!B.is_zero()) {
                if(2 * B.deg() > A.deg()) {
                    auto [a, Tr] = half_gcd(A, B);
                    concat(ak, a);
                    trs.push_back(Tr);
                    tie(A, B) = trs.back().adj().apply(A, B);
                } else {
                    auto [a, R] = A.divmod(B);
                    ak.push_back(a);
                    trs.emplace_back(a);
                    tie(A, B) = make_pair(B, R);
                }
            }
            trs.emplace_back(T(1), T(0), T(0), T(1));
            while(trs.size() >= 2) {
                trs[trs.size() - 2] = trs[trs.size() - 2] * trs[trs.size() - 1];
                trs.pop_back();
            }
            return {ak, trs.back()};
        }
                
        static poly gcd(poly A, poly B) {
            if(A.deg() < B.deg()) {
                return gcd(B, A);
            }
            auto [a, Tr] = full_gcd(A, B);
            return Tr.d * A - Tr.b * B;

        
        // Returns the characteristic polynomial
        // of the minimum linear recurrence for the sequence
        poly min_rec_slow(int d) const {
            auto R1 = mod_xk(d + 1).reverse(d + 1), R2 = xk(d + 1);
            auto Q1 = poly(T(1)), Q2 = poly(T(0));
            while(!R2.is_zero()) {
                auto [a, nR] = R1.divmod(R2); // R1 = a*R2 + nR, deg nR < deg R2
                tie(R1, R2) = make_tuple(R2, nR);
                tie(Q1, Q2) = make_tuple(Q2, Q1 + a * Q2);
                if(R2.deg() < Q2.deg()) {
                    return Q2 / Q2.lead();
                }
            }
            assert(0);
        }
        
        static transform convergent(auto L, auto R) { // computes product on [L, R)
            if(R - L == 1) {
                return transform(*L);
            } else {
                int s = 0;
                for(int i = 0; i < R - L; i++) {
                    s += L[i].a.size();
                }
                int c = 0;
                for(int i = 0; i < R - L; i++) {
                    c += L[i].a.size();
                    if(2 * c > s) {
                        return convergent(L, L + i) * convergent(L + i, R);
                    }
                }
                assert(0);
            }
        }
        
        poly min_rec(int d) const {
            if(d < magic) {
                return min_rec_slow(d);
            }
            auto R2 = mod_xk(d + 1).reverse(d + 1), R1 = xk(d + 1);
            if(R2.is_zero()) {
                return poly(1);
            }
            auto [a, Tr] = full_gcd(R1, R2);
            int dr = (d + 1) - a[0].deg();
            int dp = 0;
            for(size_t i = 0; i + 1 < a.size(); i++) {
                dr -= a[i + 1].deg();
                dp += a[i].deg();
                if(dr < dp) {
                    auto ans = convergent(begin(a), begin(a) + i + 1);
                    return ans.a / ans.a.lead();
                }
            }
            auto ans = convergent(begin(a), end(a));
            return ans.a / ans.a.lead();
        }
        
        // calculate inv to *this modulo t
        // quadratic complexity
        optional<poly> inv_mod_slow(poly const& t) const {
            auto R1 = *this, R2 = t;
            auto Q1 = poly(T(1)), Q2 = poly(T(0));
            int k = 0;
            while(!R2.is_zero()) {
                k ^= 1;
                auto [a, nR] = R1.divmod(R2);
                tie(R1, R2) = make_tuple(R2, nR);
                tie(Q1, Q2) = make_tuple(Q2, Q1 + a * Q2);
            }
            if(R1.deg() > 0) {
                return nullopt;
            } else {
                return (k ? -Q1 : Q1) / R1[0];
            }
        }
        
        optional<poly> inv_mod(poly const &t) const {
            assert(!t.is_zero());
            if(false && min(deg(), t.deg()) < magic) {
                return inv_mod_slow(t);
            }
            auto A = t, B = *this % t;
            auto [a, Tr] = full_gcd(A, B);
            auto g = Tr.d * A - Tr.b * B;
            if(g.deg() != 0) {
                return nullopt;
            }
            return -Tr.b / g[0];
        };
        
        poly operator / (const poly &t) const {return divmod(t).first;}
        poly operator % (const poly &t) const {return divmod(t).second;}
        poly operator /= (const poly &t) {return *this = divmod(t).first;}
        poly operator %= (const poly &t) {return *this = divmod(t).second;}
        poly operator *= (const T &x) {
            for(auto &it: a) {
                it *= x;
            }
            normalize();
            return *this;
        }
        poly operator /= (const T &x) {
            return *this *= x.inv();
        }
        poly operator * (const T &x) const {return poly(*this) *= x;}
        poly operator / (const T &x) const {return poly(*this) /= x;}
        
        poly conj() const { // A(x) -> A(-x)
            auto res = *this;
            for(int i = 1; i <= deg(); i += 2) {
                res.a[i] = -res[i];
            }
            return res;
        }
        
        void print(int n) const {
            for(int i = 0; i < n; i++) {
                cout << (*this)[i] << ' ';
            }
            cout << "\n";
        }
        
        void print() const {
            print(deg() + 1);
        }
        
        T eval(T x) const { // evaluates in single point x
            T res(0);
            for(int i = deg(); i >= 0; i--) {
                res *= x;
                res += a[i];
            }
            return res;
        }
        
        T lead() const { // leading coefficient
            assert(!is_zero());
            return a.back();
        }
        
        int deg() const { // degree, -1 for P(x) = 0
            return (int)a.size() - 1;
        }
        
        bool is_zero() const {
            return a.empty();
        }
        
        T operator [](int idx) const {
            return idx < 0 || idx > deg() ? T(0) : a[idx];
        }
        
        T& coef(size_t idx) { // mutable reference at coefficient
            return a[idx];
        }
        
        bool operator == (const poly &t) const {return a == t.a;}
        bool operator != (const poly &t) const {return a != t.a;}
        
        poly deriv(int k = 1) { // calculate derivative
            if(deg() + 1 < k) {
                return poly(T(0));
            }
            vector<T> res(deg() + 1 - k);
            for(int i = k; i <= deg(); i++) {
                res[i - k] = fact<T>(i) * rfact<T>(i - k) * a[i];
            }
            return res;
        }
        
        poly integr() { // calculate integral with C = 0
            vector<T> res(deg() + 2);
            for(int i = 0; i <= deg(); i++) {
                res[i + 1] = a[i] * small_inv<T>(i + 1);
            }
            return res;
        }
        
        size_t trailing_xk() const { // Let p(x) = x^k * t(x), return k
            if(is_zero()) {
                return -1;
            }
            int res = 0;
            while(a[res] == T(0)) {
                res++;
            }
            return res;
        }
        
        poly log(size_t n) { // calculate log p(x) mod x^n
            assert(a[0] == T(1));
            return (deriv().mod_xk(n) * inv(n)).integr().mod_xk(n);
        }
        
        poly exp(size_t n) { // calculate exp p(x) mod x^n
            if(is_zero()) {
                return T(1);
            }
            assert(a[0] == T(0));
            poly ans = T(1);
            size_t a = 1;
            while(a < n) {
                poly C = ans.log(2 * a).div_xk(a) - substr(a, 2 * a);
                ans -= (ans * C).mod_xk(a).mul_xk(a);
                a *= 2;
            }
            return ans.mod_xk(n);
        }
        
        poly pow_bin(int64_t k, size_t n) { // O(n log n log k)
            if(k == 0) {
                return poly(1).mod_xk(n);
            } else {
                auto t = pow(k / 2, n);
                t = (t * t).mod_xk(n);
                return (k % 2 ? *this * t : t).mod_xk(n);
            }
        }
        
        // Do not compute inverse from scratch
        poly powmod_hint(int64_t k, poly const& md, poly const& mdinv) {
            if(k == 0) {
                return poly(1);
            } else {
                auto t = powmod_hint(k / 2, md, mdinv);
                t = (t * t).divmod_hint(md, mdinv).second;
                if(k % 2) {
                    t = (t * *this).divmod_hint(md, mdinv).second;
                }
                return t;
            }
        }

        poly circular_closure(size_t m) const {
            if(deg() == -1) {
                return *this;
            }
            auto t = *this;
            for(size_t i = t.deg(); i >= m; i--) {
                t.a[i - m] += t.a[i];
            }
            t.a.resize(min(t.a.size(), m));
            return t;
        }

        static poly mul_circular(poly const& a, poly const& b, size_t m) {
            return (a.circular_closure(m) * b.circular_closure(m)).circular_closure(m);
        }

        poly powmod_circular(int64_t k, size_t m) {
            if(k == 0) {
                return poly(1);
            } else {
                auto t = powmod_circular(k / 2, m);
                t = mul_circular(t, t, m);
                if(k % 2) {
                    t = mul_circular(t, *this, m);
                }
                return t;
            }
        }
        
        poly powmod(int64_t k, poly const& md) {
            int d = md.deg();
            if(d == -1) {
                return k ? *this : poly(T(1));
            }
            if(md == xk(d)) {
                return pow(k, d);
            }
            if(md == xk(d) - poly(T(1))) {
                return powmod_circular(k, d);
            }
            auto mdinv = md.reverse().inv(md.deg() + 1);
            return powmod_hint(k, md, mdinv);
        }
        
        // O(d * n) with the derivative trick from
        // https://codeforces.com/blog/entry/73947?#comment-581173
        poly pow_dn(int64_t k, size_t n) {
            if(n == 0) {
                return poly(T(0));
            }
            assert((*this)[0] != T(0));
            vector<T> Q(n);
            Q[0] = bpow(a[0], k);
            auto a0inv = a[0].inv();
            for(int i = 1; i < (int)n; i++) {
                for(int j = 1; j <= min(deg(), i); j++) {
                    Q[i] += a[j] * Q[i - j] * (T(k) * T(j) - T(i - j));
                }
                Q[i] *= small_inv<T>(i) * a0inv;
            }
            return Q;
        }
        
        // calculate p^k(n) mod x^n in O(n log n)
        // might be quite slow due to high constant
        poly pow(int64_t k, size_t n) {
            if(is_zero()) {
                return k ? *this : poly(1);
            }
            int i = trailing_xk();
            if(i > 0) {
                return k >= int64_t(n + i - 1) / i ? poly(T(0)) : div_xk(i).pow(k, n - i * k).mul_xk(i * k);
            }
            if(min(deg(), (int)n) <= magic) {
                return pow_dn(k, n);
            }
            if(k <= magic) {
                return pow_bin(k, n);
            }
            T j = a[i];
            poly t = *this / j;
            return bpow(j, k) * (t.log(n) * T(k)).exp(n).mod_xk(n);
        }
        
        // returns nullopt if undefined
        optional<poly> sqrt(size_t n) const {
            if(is_zero()) {
                return *this;
            }
            int i = trailing_xk();
            if(i % 2) {
                return nullopt;
            } else if(i > 0) {
                auto ans = div_xk(i).sqrt(n - i / 2);
                return ans ? ans->mul_xk(i / 2) : ans;
            }
            auto st = (*this)[0].sqrt();
            if(st) {
                poly ans = *st;
                size_t a = 1;
                while(a < n) {
                    a *= 2;
                    ans -= (ans - mod_xk(a) * ans.inv(a)).mod_xk(a) / 2;
                }
                return ans.mod_xk(n);
            }
            return nullopt;
        }
        
        poly mulx(T a) const { // component-wise multiplication with a^k
            T cur = 1;
            poly res(*this);
            for(int i = 0; i <= deg(); i++) {
                res.coef(i) *= cur;
                cur *= a;
            }
            return res;
        }

        poly mulx_sq(T a) const { // component-wise multiplication with a^{k choose 2}
            T cur = 1, total = 1;
            poly res(*this);
            for(int i = 0; i <= deg(); i++) {
                res.coef(i) *= total;
                cur *= a;
                total *= cur;
            }
            return res;
        }

        // be mindful of maxn, as the function
        // requires multiplying polynomials of size deg() and n+deg()!
        poly chirpz(T z, int n) const { // P(1), P(z), P(z^2), ..., P(z^(n-1))
            if(is_zero()) {
                return vector<T>(n, 0);
            }
            if(z == T(0)) {
                vector<T> ans(n, (*this)[0]);
                if(n > 0) {
                    ans[0] = accumulate(begin(a), end(a), T(0));
                }
                return ans;
            }
            auto A = mulx_sq(z.inv());
            auto B = ones(n+deg()).mulx_sq(z);
            return semicorr(B, A).mod_xk(n).mulx_sq(z.inv());
        }

        // res[i] = prod_{1 <= j <= i} 1/(1 - z^j)
        static auto _1mzk_prod_inv(T z, int n) {
            vector<T> res(n, 1), zk(n);
            zk[0] = 1;
            for(int i = 1; i < n; i++) {
                zk[i] = zk[i - 1] * z;
                res[i] = res[i - 1] * (T(1) - zk[i]);
            }
            res.back() = res.back().inv();
            for(int i = n - 2; i >= 0; i--) {
                res[i] = (T(1) - zk[i+1]) * res[i+1];
            }
            return res;
        }
        
        // prod_{0 <= j < n} (1 - z^j x)
        static auto _1mzkx_prod(T z, int n) {
            if(n == 1) {
                return poly(vector<T>{1, -1});
            } else {
                auto t = _1mzkx_prod(z, n / 2);
                t *= t.mulx(bpow(z, n / 2));
                if(n % 2) {
                    t *= poly(vector<T>{1, -bpow(z, n - 1)});
                }
                return t;
            }
        }

        poly chirpz_inverse(T z, int n) const { // P(1), P(z), P(z^2), ..., P(z^(n-1))
            if(is_zero()) {
                return {};
            }
            if(z == T(0)) {
                if(n == 1) {
                    return *this;
                } else {
                    return vector{(*this)[1], (*this)[0] - (*this)[1]};
                }
            }
            vector<T> y(n);
            for(int i = 0; i < n; i++) {
                y[i] = (*this)[i];
            }
            auto prods_pos = _1mzk_prod_inv(z, n);
            auto prods_neg = _1mzk_prod_inv(z.inv(), n);

            T zn = bpow(z, n-1).inv();
            T znk = 1;
            for(int i = 0; i < n; i++) {
                y[i] *= znk * prods_neg[i] * prods_pos[(n - 1) - i];
                znk *= zn;
            }

            poly p_over_q = poly(y).chirpz(z, n);
            poly q = _1mzkx_prod(z, n);

            return (p_over_q * q).mod_xk(n).reverse(n);
        }

        static poly build(vector<poly> &res, int v, auto L, auto R) { // builds evaluation tree for (x-a1)(x-a2)...(x-an)
            if(R - L == 1) {
                return res[v] = vector<T>{-*L, 1};
            } else {
                auto M = L + (R - L) / 2;
                return res[v] = build(res, 2 * v, L, M) * build(res, 2 * v + 1, M, R);
            }
        }

        poly to_newton(vector<poly> &tree, int v, auto l, auto r) {
            if(r - l == 1) {
                return *this;
            } else {
                auto m = l + (r - l) / 2;
                auto A = (*this % tree[2 * v]).to_newton(tree, 2 * v, l, m);
                auto B = (*this / tree[2 * v]).to_newton(tree, 2 * v + 1, m, r);
                return A + B.mul_xk(m - l);
            }
        }

        poly to_newton(vector<T> p) {
            if(is_zero()) {
                return *this;
            }
            int n = p.size();
            vector<poly> tree(4 * n);
            build(tree, 1, begin(p), end(p));
            return to_newton(tree, 1, begin(p), end(p));
        }

        vector<T> eval(vector<poly> &tree, int v, auto l, auto r) { // auxiliary evaluation function
            if(r - l == 1) {
                return {eval(*l)};
            } else {
                auto m = l + (r - l) / 2;
                auto A = (*this % tree[2 * v]).eval(tree, 2 * v, l, m);
                auto B = (*this % tree[2 * v + 1]).eval(tree, 2 * v + 1, m, r);
                A.insert(end(A), begin(B), end(B));
                return A;
            }
        }
        
        vector<T> eval(vector<T> x) { // evaluate polynomial in (x1, ..., xn)
            int n = x.size();
            if(is_zero()) {
                return vector<T>(n, T(0));
            }
            vector<poly> tree(4 * n);
            build(tree, 1, begin(x), end(x));
            return eval(tree, 1, begin(x), end(x));
        }
        
        poly inter(vector<poly> &tree, int v, auto ly, auto ry) { // auxiliary interpolation function
            if(ry - ly == 1) {
                return {*ly / a[0]};
            } else {
                auto my = ly + (ry - ly) / 2;
                auto A = (*this % tree[2 * v]).inter(tree, 2 * v, ly, my);
                auto B = (*this % tree[2 * v + 1]).inter(tree, 2 * v + 1, my, ry);
                return A * tree[2 * v + 1] + B * tree[2 * v];
            }
        }
        
        static auto inter(vector<T> x, vector<T> y) { // interpolates minimum polynomial from (xi, yi) pairs
            int n = x.size();
            vector<poly> tree(4 * n);
            return build(tree, 1, begin(x), end(x)).deriv().inter(tree, 1, begin(y), end(y));
        }

        static auto resultant(poly a, poly b) { // computes resultant of a and b
            if(b.is_zero()) {
                return 0;
            } else if(b.deg() == 0) {
                return bpow(b.lead(), a.deg());
            } else {
                int pw = a.deg();
                a %= b;
                pw -= a.deg();
                auto mul = bpow(b.lead(), pw) * T((b.deg() & a.deg() & 1) ? -1 : 1);
                auto ans = resultant(b, a);
                return ans * mul;
            }
        }
                
        static poly xk(size_t n) { // P(x) = x^n
            return poly(T(1)).mul_xk(n);
        }
        
        static poly ones(size_t n) { // P(x) = 1 + x + ... + x^{n-1} 
            return vector<T>(n, 1);
        }
        
        static poly expx(size_t n) { // P(x) = e^x (mod x^n)
            return ones(n).borel();
        }

        static poly log1px(size_t n) { // P(x) = log(1+x) (mod x^n)
            vector<T> coeffs(n, 0);
            for(size_t i = 1; i < n; i++) {
                coeffs[i] = (i & 1 ? T(i).inv() : -T(i).inv());
            }
            return coeffs;
        }

        static poly log1mx(size_t n) { // P(x) = log(1-x) (mod x^n)
            return -ones(n).integr();
        }
        
        // [x^k] (a corr b) = sum_{i} a{(k-m)+i}*bi
        static poly corr(poly a, poly b) { // cross-correlation
            return a * b.reverse();
        }

        // [x^k] (a semicorr b) = sum_i a{i+k} * b{i}
        static poly semicorr(poly a, poly b) {
            return corr(a, b).div_xk(b.deg());
        }
        
        poly invborel() const { // ak *= k!
            auto res = *this;
            for(int i = 0; i <= deg(); i++) {
                res.coef(i) *= fact<T>(i);
            }
            return res;
        }
        
        poly borel() const { // ak /= k!
            auto res = *this;
            for(int i = 0; i <= deg(); i++) {
                res.coef(i) *= rfact<T>(i);
            }
            return res;
        }
        
        poly shift(T a) const { // P(x + a)
            return semicorr(invborel(), expx(deg() + 1).mulx(a)).borel();
        }
        
        poly x2() { // P(x) -> P(x^2)
            vector<T> res(2 * a.size());
            for(size_t i = 0; i < a.size(); i++) {
                res[2 * i] = a[i];
            }
            return res;
        }
        
        // Return {P0, P1}, where P(x) = P0(x) + xP1(x)
        pair<poly, poly> bisect() const {
            vector<T> res[2];
            res[0].reserve(deg() / 2 + 1);
            res[1].reserve(deg() / 2 + 1);
            for(int i = 0; i <= deg(); i++) {
                res[i % 2].push_back(a[i]);
            }
            return {res[0], res[1]};
        }
        
        // Find [x^k] P / Q
        static T kth_rec(poly P, poly Q, int64_t k) {
            while(k > Q.deg()) {
                int n = Q.a.size();
                auto [Q0, Q1] = Q.mulx(-1).bisect();
                auto [P0, P1] = P.bisect();
                
                int N = fft::com_size((n + 1) / 2, (n + 1) / 2);
                
                auto Q0f = fft::dft(Q0.a, N);
                auto Q1f = fft::dft(Q1.a, N);
                auto P0f = fft::dft(P0.a, N);
                auto P1f = fft::dft(P1.a, N);
                
                if(k % 2) {
                    P = poly(Q0f * P1f) + poly(Q1f * P0f);
                } else {
                    P = poly(Q0f * P0f) + poly(Q1f * P1f).mul_xk(1);
                }
                Q = poly(Q0f * Q0f) - poly(Q1f * Q1f).mul_xk(1);
                k /= 2;
            }
            return (P * Q.inv(Q.deg() + 1))[k];
        }
        
        poly inv(int n) const { // get inverse series mod x^n
            auto Q = mod_xk(n);
            if(n == 1) {
                return Q[0].inv();
            }
            // Q(-x) = P0(x^2) + xP1(x^2)
            auto [P0, P1] = Q.mulx(-1).bisect();
            
            int N = fft::com_size((n + 1) / 2, (n + 1) / 2);
            
            auto P0f = fft::dft(P0.a, N);
            auto P1f = fft::dft(P1.a, N);
            
            auto TTf = fft::dft(( // Q(x)*Q(-x) = Q0(x^2)^2 - x^2 Q1(x^2)^2
                poly(P0f * P0f) - poly(P1f * P1f).mul_xk(1)
            ).inv((n + 1) / 2).a, N);
            
            return (
                poly(P0f * TTf).x2() + poly(P1f * TTf).x2().mul_xk(1)
            ).mod_xk(n);
        }
        
        // compute A(B(x)) mod x^n in O(n^2)
        static poly compose(poly A, poly B, int n) {
            int q = std::sqrt(n);
            vector<poly> Bk(q);
            auto Bq = B.pow(q, n);
            Bk[0] = poly(T(1));
            for(int i = 1; i < q; i++) {
                Bk[i] = (Bk[i - 1] * B).mod_xk(n);
            }
            poly Bqk(1);
            poly ans;
            for(int i = 0; i <= n / q; i++) {
                poly cur;
                for(int j = 0; j < q; j++) {
                    cur += Bk[j] * A[i * q + j];
                }
                ans += (Bqk * cur).mod_xk(n);
                Bqk = (Bqk * Bq).mod_xk(n);
            }
            return ans;
        }
        
        // compute A(B(x)) mod x^n in O(sqrt(pqn log^3 n))
        // preferrable when p = deg A and q = deg B
        // are much less than n
        static poly compose_large(poly A, poly B, int n) {
            if(B[0] != T(0)) {
                return compose_large(A.shift(B[0]), B - B[0], n);
            }
            
            int q = std::sqrt(n);
            auto [B0, B1] = make_pair(B.mod_xk(q), B.div_xk(q));
            
            B0 = B0.div_xk(1);
            vector<poly> pw(A.deg() + 1);
            auto getpow = [&](int k) {
                return pw[k].is_zero() ? pw[k] = B0.pow(k, n - k) : pw[k];
            };
            
            function<poly(poly const&, int, int)> compose_dac = [&getpow, &compose_dac](poly const& f, int m, int N) {
                if(f.deg() <= 0) {
                    return f;
                }
                int k = m / 2;
                auto [f0, f1] = make_pair(f.mod_xk(k), f.div_xk(k));
                auto [A, B] = make_pair(compose_dac(f0, k, N), compose_dac(f1, m - k, N - k));
                return (A + (B.mod_xk(N - k) * getpow(k).mod_xk(N - k)).mul_xk(k)).mod_xk(N);
            };
            
            int r = n / q;
            auto Ar = A.deriv(r);
            auto AB0 = compose_dac(Ar, Ar.deg() + 1, n);
            
            auto Bd = B0.mul_xk(1).deriv();
            
            poly ans = T(0);
            
            vector<poly> B1p(r + 1);
            B1p[0] = poly(T(1));
            for(int i = 1; i <= r; i++) {
                B1p[i] = (B1p[i - 1] * B1.mod_xk(n - i * q)).mod_xk(n - i * q);
            }
            while(r >= 0) {
                ans += (AB0.mod_xk(n - r * q) * rfact<T>(r) * B1p[r]).mul_xk(r * q).mod_xk(n);
                r--;
                if(r >= 0) {
                    AB0 = ((AB0 * Bd).integr() + A[r] * fact<T>(r)).mod_xk(n);
                }
            }
            
            return ans;
        }
    };
    
    static auto operator * (const auto& a, const poly<auto>& b) {
        return b * a;
    }
};
