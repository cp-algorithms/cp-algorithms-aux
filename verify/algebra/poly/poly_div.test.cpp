// @brief Division of Polynomials
#define PROBLEM "https://judge.yosupo.jp/problem/division_of_polynomials"
#include "cp-algo/algebra/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;

const int mod = 998244353;
using base = modular<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, m;
    cin >> n >> m;
    vector<base> a(n), b(m);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    copy_n(istream_iterator<base>(cin), m, begin(b));
    auto [q, r] = polyn(a).divmod(b);
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
