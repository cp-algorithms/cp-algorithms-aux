// @brief Product of Polynomial Sequence
#define PROBLEM "https://judge.yosupo.jp/problem/product_of_polynomial_sequence"
#pragma GCC optimize("O3,unroll-loops")
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int N;
    cin >> N;
    vector<polyn> polys(N);
    multiset<polyn, decltype([](polyn const& a, polyn const& b){
        return a.deg() < b.deg();
    })> que = {polyn(1)};
    int D = 0;
    for(int i = 0; i < N; i++) {
        int d;
        cin >> d;
        D += d;
        polyn::Vector a(d + 1);
        copy_n(istream_iterator<base>(cin), d + 1, begin(a));
        que.insert(polyn(a));
    }
    while(que.size() > 1) {
        auto A = *begin(que);
        que.erase(begin(que));
        auto B = *begin(que);
        que.erase(begin(que));
        que.insert(A * B);
    }
    begin(que)->print(D + 1);
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
