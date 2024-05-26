// @brief Wildcard Pattern Matching
#define PROBLEM "https://judge.yosupo.jp/problem/wildcard_pattern_matching"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

using base = complex<double>;
using polyn = poly_t<base>;

string matches(string const& A, string const& B, char wild = '*') {
    static base c_to_int[2][26];
    static bool init = false;
    if(!init) {
        init = true;
        for(int i = 0; i < 26; i++) {
            c_to_int[0][i] = polar(1., (double)cp_algo::random::rng());
            c_to_int[1][i] = conj(c_to_int[0][i]);
        }
    }
    string ST[2] = {A, B};
    polyn P[2];
    for(int i: {0, 1}) {
        vector<base> coeffs(size(ST[i]));
        for(size_t k = 0; k < size(ST[i]); k++) {
            coeffs[k] = base(ST[i][k] == wild ? 0 : c_to_int[i][ST[i][k] - 'a']);
        }
        P[i] = coeffs;
    }
    auto dist0 = polyn::inner_semicorr(P[0], P[1]);
    string ans(size(ST[0]) - size(ST[1]) + 1, '0');
    for(size_t j = 0; j <= size(ans); j++) {
        ans[j] = '0' + (
            abs(dist0[j].imag()) < 1e-8 && abs(dist0[j].real() - round(dist0[j].real())) < 1e-8
        );
    }
    return ans;
}

void solve() {
    string a, b;
    cin >> a >> b;
    cout << matches(a, b) << "\n";
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
