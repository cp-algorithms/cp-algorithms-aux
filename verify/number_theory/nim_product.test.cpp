// @brief Nim Product
#define PROBLEM "https://judge.yosupo.jp/problem/nim_product_64"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt,pclmul")
#define CP_ALGO_CHECKPOINT
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include <immintrin.h>
#include "cp-algo/util/big_alloc.hpp"
#include "cp-algo/util/checkpoint.hpp"
#include "cp-algo/number_theory/nimber.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math::nimber;

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t;
    cin >> t;
    cp_algo::big_vector<uint64_t> A(t), B(t);
    for (int i = 0; i < t; i++) {
        cin >> A[i] >> B[i];
    }
    cp_algo::checkpoint("read");
    for (int i = 0; i < t; i++) {
        A[i] = nim_to_poly(A[i]);
        B[i] = nim_to_poly(B[i]);
    }
    cp_algo::checkpoint("to_poly");
    for (int i = 0; i < t; i++) {
        A[i] = reduce_mod(clmul(A[i], B[i]));
    }
    cp_algo::checkpoint("clmul+reduce");
    for (int i = 0; i < t; i++) {
        A[i] = poly_to_nim(A[i]);
    }
    cp_algo::checkpoint("to_nim");
    for (int i = 0; i < t; i++) {
        cout << A[i] << "\n";
    }
    cp_algo::checkpoint("print");
    cp_algo::checkpoint<1>();
}
