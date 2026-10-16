#include "cp-algo/linalg/matrix.hpp"
#include "cp-algo/math/multivar.hpp"
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <sstream>

using namespace cp_algo;

template<class T> concept has_reduced_reader = requires(T value) { value.read_reduced(); };
static_assert(!has_reduced_reader<linalg::vec<double>>);
static_assert(!has_reduced_reader<linalg::matrix<double, linalg::vec<double>>>);

template<class T> void check() {
    constexpr int64_t m = 998244353;
    std::istringstream input("-1 998244354 0 998244352 1 0 2 3 998244352 1 0 2");
    auto *original = std::cin.rdbuf(input.rdbuf());
    linalg::vec<T> normal(size_t(2));
    normal.read();
    assert(normal[0].getr() == m - 1 && normal[1].getr() == 1);
    linalg::vec<T> vector(size_t(2));
    vector.read_reduced();
    assert(vector[0].getr() == 0 && vector[1].getr() == m - 1);
    using Row = std::conditional_t<(T::bits >= 64), linalg::modint_vec<T>, linalg::vec<T>>;
    linalg::matrix<T, Row> matrix(2, 2);
    matrix.read_reduced();
    assert(matrix[0][0].getr() == 1 && matrix[0][1].getr() == 0);
    assert(matrix[1][0].getr() == 2 && matrix[1][1].getr() == 3);
    math::fft::multivar<T> polynomial(std::array<size_t, 2>{2, 2});
    polynomial.read_reduced();
    constexpr std::array<int64_t, 4> expected{m - 1, 1, 0, 2};
    for(size_t i = 0; i < 4; i++) {
        assert(polynomial.data[i].getr() == expected[i]);
    }
    assert(!std::cin.fail());
    std::cin.rdbuf(original);

    std::istringstream failed("2 invalid");
    original = std::cin.rdbuf(failed.rdbuf());
    vector[0] = 5;
    vector[1] = 7;
    vector.read_reduced();
    assert(std::cin.fail());
    assert(vector[0].getr() == 2 && vector[1].getr() == 7);
    std::cin.rdbuf(original);
    std::cin.clear();
}

int main() {
    check<math::modint<998244353>>();
    check<math::modint<998244353LL>>();
    math::dynamic_modint<>::with_mod(998244353, [] { check<math::dynamic_modint<>>(); });
    math::dynamic_modint<int64_t>::with_mod(998244353, [] { check<math::dynamic_modint<int64_t>>(); });
}
