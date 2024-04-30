// @brief Characteristic Polynomial
#define PROBLEM "https://judge.yosupo.jp/problem/characteristic_polynomial"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/algebra/poly.hpp"
#include "cp-algo/linalg/matrix.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;
using namespace cp_algo::linalg;

const int mod = 998244353;
using base = modular<mod>;
using polyn = poly_t<base>;

template<bool reduce = false>
auto frobenius_basis(auto const& A) {
    assert(A.n() == A.m());
    size_t n = A.n();
    struct krylov {
        vector<vec<base>> basis; 
        poly_t<base> rec;
    };
    vector<krylov> blocks;
    vector<vec<base>> reduced;
    while(size(reduced) < n) {
        auto generate_block = [&](auto x) {
            krylov block;
            while(true) {
                vec<base> y = x | vec<base>::ei(n + 1, size(reduced));
                for(auto &it: reduced) {
                    y.reduce_by(it);
                }
                y.normalize();
                if(vec<base>(y[slice(0, n, 1)]) == vec<base>(n)) {
                    block.rec = vector<base>(
                        begin(y) + n + size(reduced) - size(block.basis),
                        begin(y) + n + size(reduced) + 1
                    );
                    return std::pair{block, vec<base>(y[slice(n, n, 1)])};
                } else {
                    block.basis.push_back(x);
                    reduced.push_back(y);
                    x = A.apply(x);
                }
            }
        };
        auto [block, full_rec] = generate_block(vec<base>::random(n));
        if constexpr (reduce) {
            if(vec<base>(full_rec[slice(0, size(reduced), 1)]) != vec<base>(size(reduced))) {
                auto x = block.basis[0];
                size_t start = 0;
                for(auto &[basis, rec]: blocks) {
                    polyn cur_rec = vector<base>(
                        begin(full_rec) + start, begin(full_rec) + start + rec.deg()
                    );
                    auto shift = cur_rec / block.rec;
                    for(int j = 0; j <= shift.deg(); j++) {
                        x.add_scaled(basis[j], shift[j]);
                    }
                    start += rec.deg();
                }
                reduced.erase(begin(reduced) + start, end(reduced));
                tie(block, full_rec) = generate_block(x.normalize());
            }
        }
        blocks.push_back(block);
    }
    return blocks;
}

void solve() {
    size_t n;
    cin >> n;
    matrix<base> A(n);
    A.read();
    auto blocks = frobenius_basis(A);
    polyn res(1);
    for(auto &[basis, rec]: blocks) {
        res *= rec;
    }
    res.print();
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