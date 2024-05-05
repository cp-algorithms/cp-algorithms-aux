// @brief Polynomial Root Finding
#define PROBLEM "https://judge.yosupo.jp/problem/polynomial_root_finding"
#include "cp-algo/algebra/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;
using namespace cp_algo::random;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void find_roots_impl(polyn const& p, vector<base> &res) {
    if(p.deg() == 1) {
        res.push_back(-p[0] / p[1]);
    } else if(p.deg() > 1) {
        auto A = polyn::gcd(polyn({rng(), 1}).powmod((mod - 1) / 2, p) - base(1), polyn(p));
        find_roots_impl(A, res);
        find_roots_impl(p / A, res);
    }
}

auto find_roots(polyn const& p) {
    vector<base> res;
    if(p[0] == 0) {
        res.push_back(0);
    }
    auto g = polyn::xk(1).powmod(mod - 1, p);
    find_roots_impl(polyn::gcd(g - base(1), polyn(p)), res);
    return res;
}

void solve() {
    int n;
    cin >> n;
    vector<base> f(n+1);
    copy_n(istream_iterator<base>(cin), n+1, begin(f));
    polyn res = find_roots(f);
    cout << res.deg() + 1 << "\n";
    res.print();
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