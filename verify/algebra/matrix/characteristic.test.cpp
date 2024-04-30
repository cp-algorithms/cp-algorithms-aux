// @brief Characteristic Polynomial
#define PROBLEM "https://judge.yosupo.jp/problem/characteristic_polynomial"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/algebra/poly.hpp"
#include "cp-algo/linalg/matrix.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;
using namespace cp_algo::linalg;

const int mod = 998244353;
using base = modular<mod>;
using polyn = poly_t<base>;

void solve() {
    size_t n;
    cin >> n;
    matrix<base> A(n);
    A.read();
    vector<vec<base>> basis;
    auto x = vec<base>::random(n);
    size_t degree = 0;
    polyn ans = base(1);
    while(size(basis) <= n) {
        auto y = x | vec<base>::ei(n + 1, size(basis));
        for(auto &it: basis) {
            y.reduce_by(it);
        }
        y.normalize();
        if(vec<base>(y[slice(0, n, 1)]) == vec<base>(n)) {
            vector<base> cur(begin(y) + n + size(basis) - degree,
                             begin(y) + n + size(basis) + 1);
            ans *= polyn(cur);
            degree = 0;
            if(size(basis) < n) {
                x = vec<base>::random(n);
            } else {
                break;
            }
        } else {
            basis.push_back(y);
            x = A.apply(x);
            degree++;
        }
    }
    ans.print();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    //cin >> t;
    while(t--) {
        solve();
    }
}
