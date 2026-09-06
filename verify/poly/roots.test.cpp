// @brief Polynomial Root Finding
#define PROBLEM "https://judge.yosupo.jp/problem/polynomial_root_finding"

#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/euclid.hpp"
#include "cp-algo/math/poly/powmod.hpp"

using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::random;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void find_roots_impl(polyn const& p, polyn::Vector &res) {
    if(p.deg() == 1) {
        res.push_back(-p[0] / p[1]);
    } else if(p.deg() > 1) {
        auto A = gcd(powmod(polyn(polyn::Vector{(base)rng(), 1}), (mod - 1) / 2, p) - base(1), polyn(p));
        find_roots_impl(A, res);
        find_roots_impl(p / A, res);
    }
}

auto find_roots(polyn const& p) {
    polyn::Vector res;
    if(p[0] == 0) {
        res.push_back(0);
    }
    auto g = powmod(polyn::xk(1), mod - 1, p);
    find_roots_impl(gcd(g - base(1), polyn(p)), res);
    return res;
}

void solve() {
    int n;
    cin >> n;
    polyn::Vector f(n+1);
    for(auto &it: f) {cin >> it;}
    auto res = find_roots(f);
    cout << res.size() << "\n";
    for(auto &it: res) {cout << it << ' ';}
    cout << "\n";
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
