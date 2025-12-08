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

simd_target auto is_integer(auto a) {
    static const ftype eps = 1e-9;
    return cp_algo::abs(a - cp_algo::round(a)) < eps;
}

simd_target string matches(string const& A, string const& B, char wild = '*') {
    static ftype project[2][128];
    static bool init = false;
    if(!init) {
        init = true;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(.5, 2.);
        for(int i = 0; i < 128; i++) {
            ftype x = dis(gen);
            project[0][i] = x;
            project[1][i] = 1. / x;
        }
    }
    project[0][(int)wild] = project[1][(int)wild] = 0;
    vector<cvector> P;
    P.emplace_back((size(A) + 1) / 2);
    P.emplace_back((size(A) + 1) / 2);
    auto N = P[0].size();
    auto assign = [&](int z) {
        return [&, z](auto ic) {
            auto [i, c] = ic;
            if(i < (int)N) {
                real(P[z].r[i / fft::flen])[i % fft::flen] = project[z][(int)c];
            } else {
                i -= N;
                imag(P[z].r[i / fft::flen])[i % fft::flen] = project[z][(int)c];
            }
        };
    };
    ranges::for_each(A | views::enumerate, assign(0));
    ranges::for_each(B | views::reverse | views::enumerate, assign(1));
    cp_algo::checkpoint("cvector fill");
    semicorr(P[0], P[1]);
    string ans(2 * size(P[0]), '0');
    auto start = (ssize(B) - 1) / fft::flen * fft::flen;
    for(auto j = start; j < size(ans); j += fft::flen) {
        decltype(is_integer(real(P[0].at(j)))) check;
        if(j < N) {
            check = is_integer(real(P[0].at(j)));
        } else {
            check = is_integer(imag(P[0].at(j - N)));
        }
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
