// @brief Represent A Number As Two Square Sum
#define PROBLEM "https://judge.yosupo.jp/problem/two_square_sum"
#pragma GCC optimize("O3,unroll-loops")
#include "cp-algo/number_theory/two_squares.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

void solve() {
	int64_t n;
	cin >> n;
	auto res = two_squares_all(n);
	cout << size(res) << "\n";
	for(auto p: res) {
		cout << p.real() << ' ' << p.imag() << "\n";
	}
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    cin >> t;
	while(t--) {
		solve();
	}
}
