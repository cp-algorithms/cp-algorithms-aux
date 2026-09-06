// @brief Product of Polynomial Sequence
#define PROBLEM "https://judge.yosupo.jp/problem/product_of_polynomial_sequence"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/base.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int N;
    cin >> N;
    vector<optional<polyn>> prod;
    int D = 0;
    base scalar = 1;
    for(int i = 0; i < N; i++) {
        int d;
        cin >> d;
        D += d;
        polyn::Vector a(d + 1);
        for(auto &it: a) {cin >> it;}
        polyn p(std::move(a));
        if(p.deg() == 0) {scalar *= p[0]; continue;}
        // Merge products of comparable sizes as they arrive.
        while(true) {
            auto k = std::bit_width(size_t(std::max(0, p.deg())));
            if(k >= prod.size()) {prod.resize(k + 1);}
            if(!prod[k]) {prod[k] = std::move(p); break;}
            p *= *prod[k];
            prod[k].reset();
        }
    }
    polyn ans(1);
    for(auto &p: prod) {
        if(p) {ans *= *p;}
    }
    if(scalar != base(1)) {ans *= scalar;}
    ans.print(D + 1);
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
