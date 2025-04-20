// @brief Wildcard Pattern Matching
#define PROBLEM "https://judge.yosupo.jp/problem/wildcard_pattern_matching"
#pragma GCC optimize("Ofast,unroll-loops")
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
    array ST = {&A, &B};
    vector<cvector> P;
    P.emplace_back(size(A));
    P.emplace_back(size(A));
    for(int i: {0, 1}) {
        size_t N = ST[i]->size();
        for(size_t k = 0; k < N; k++) {
            char c = ST[i]->at(k);
            size_t idx = i ? N - k - 1 : k;
            point val = c == wild ? 0 : project[i][c - 'a'];
            P[i].set(idx, val);
        }
    }
    semicorr(P[0], P[1]);
    string ans(size(A) - size(B) + 1, '0');
    for(size_t j = 0; j < size(ans); j++) {
        ans[j] = '0' + is_integer(P[0].get(size(B) - 1 + j));
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
