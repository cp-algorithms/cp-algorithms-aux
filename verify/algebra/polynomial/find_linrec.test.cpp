// @brief Find Linear Recurrence
#define PROBLEM "https://judge.yosupo.jp/problem/find_linear_recurrence"
#include "cp-algo/algebra/polynomial.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;

const int mod = 998244353;
typedef modular<mod> base;
typedef poly<base> polyn;

void solve() {
    int n;
    cin >> n;
    vector<base> a(n);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    auto Q = polyn(a).min_rec(n - 1);
    int d = Q.deg();
    cout << d << endl;
    (-Q.reverse().div_xk(1)).print(d);
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    while(t--) {
        solve();
    }
}
