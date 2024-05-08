// @brief Wildcard Pattern Matching
#define PROBLEM "https://judge.yosupo.jp/problem/wildcard_pattern_matching"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/algebra/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;

const int mod = 1e9 + 9;

using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    string ST[2];
    cin >> ST[0] >> ST[1];
    polyn P[2][3];
    for(int i: {0, 1}) {
        for(int j: {1, 2, 3}) {
            vector<base> coeffs(size(ST[i]));
            for(size_t k = 0; k < size(ST[i]); k++) {
                coeffs[k] = bpow(base(ST[i][k] == '*' ? 0 : ST[i][k] - 'a' + 1), j);
            }
            P[i][j-1] = coeffs;
        }
    }
    auto dist = polyn::semicorr(P[0][0], P[1][2])
        + polyn::semicorr(P[0][2], P[1][0])
        - 2*polyn::semicorr(P[0][1], P[1][1]);
    string ans(size(ST[0]) - size(ST[1]) + 1, '0');
    for(size_t j = 0; j <= size(ans); j++) {
        ans[j] = '0' + int(dist[j] == 0);
    }
    cout << ans << "\n";
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
