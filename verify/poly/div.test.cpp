// @brief Division of Polynomials
#define PROBLEM "https://judge.yosupo.jp/problem/division_of_polynomials"
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/div.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, m;
    cin >> n >> m;
    polyn::Vector a(n), b(m);
    for(auto &it: a) {cin >> it;}
    for(auto &it: b) {cin >> it;}
    auto [q, r] = divmod(polyn(std::move(a)), polyn(std::move(b)));
    cout << q.deg() + 1 << ' ' << r.deg() + 1 << "\n";
    q.print();
    r.print();
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
