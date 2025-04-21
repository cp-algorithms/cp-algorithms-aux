// @brief Wildcard Pattern Matching
#define PROBLEM "https://judge.yosupo.jp/problem/wildcard_pattern_matching"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include "cp-algo/math/cvector.hpp"
#include "cp-algo/random/rng.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

using fft::ftype;
using fft::point;
using fft::cvector;

void semicorr(auto &a, auto &b) {
    a.fft();
    b.fft();
    a.dot(b);
    a.ifft();
}

auto is_integer = [](point a) {
    static const double eps = 1e-8;
    return abs(imag(a)) < eps
        && abs(real(a) - round(real(a))) < eps;
};

string matches(string const& A, string const& B, char wild = '*') {
    static const int sigma = 26;
    static point project[2][sigma];
    static bool init = false;
    if(!init) {
        init = true;
        for(int i = 0; i < sigma; i++) {
            project[0][i] = cp_algo::polar(1., (ftype)cp_algo::random::rng());
            project[1][i] = conj(project[0][i]);
        }
    }
    vector<cvector> P;
    P.emplace_back(size(A));
    P.emplace_back(size(A));
    for(auto [i, c]: A | views::enumerate) {
        P[0].set(i, (c != wild) * project[0][c - 'a']);
    }
    for(auto [i, c]: B | views::reverse | views::enumerate) {
        P[1].set(i, (c != wild) * project[1][c - 'a']);
    }
    cp_algo::checkpoint("cvector fill");
    semicorr(P[0], P[1]);
    string ans(size(A) - size(B) + 1, '0');
    for(size_t j = 0; j < size(ans); j++) {
        ans[j] = '0' + is_integer(P[0].get(size(B) - 1 + j));
    }
    cp_algo::checkpoint("fill answer");
    return ans;
}

void solve() {
    string a, b;
    cin >> a >> b;
    cp_algo::checkpoint("input");
    cout << matches(a, b) << "\n";
    cp_algo::checkpoint("output");
    cp_algo::checkpoint<true>("done");
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
