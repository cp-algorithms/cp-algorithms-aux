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
using fft::vftype;
using fft::cvector;

void semicorr(auto &a, auto &b) {
    a.fft();
    b.fft();
    a.dot(b);
    a.ifft();
}

auto is_integer(auto a) {
    static const double eps = 1e-8;
    return cp_algo::abs(imag(a)) < eps
        && cp_algo::abs(real(a) - cp_algo::round(real(a))) < eps;
}

string matches(string const& A, string const& B, char wild = '*') {
    static point project[2][128];
    static bool init = false;
    if(!init) {
        init = true;
        for(int i = 0; i < 128; i++) {
            project[0][i] = cp_algo::polar(1., (ftype)cp_algo::random::rng());
            project[1][i] = conj(project[0][i]);
        }
    }
    project[0][(int)wild] = project[1][(int)wild] = 0;
    vector<cvector> P;
    P.emplace_back(size(A));
    P.emplace_back(size(A));
    for(auto [i, c]: A | views::enumerate) {
        P[0].set(i, project[0][(int)c]);
    }
    for(auto [i, c]: B | views::reverse | views::enumerate) {
        P[1].set(i, project[1][(int)c]);
    }
    cp_algo::checkpoint("cvector fill");
    semicorr(P[0], P[1]);
    string ans(size(P[0]), '0');
    auto start = (size(B) - 1) / fft::flen * fft::flen;
    for(size_t j = start; j < size(ans); j += fft::flen) {
        auto r = P[0].at(j);
        auto check = is_integer(r);
        for(int z = 0; z < 4; z++) {
            ans[j + z] ^= (bool)check[z];
        }
    }
    cp_algo::checkpoint("fill answer");
    return ans.substr(size(B) - 1, size(A) - size(B) + 1);
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
