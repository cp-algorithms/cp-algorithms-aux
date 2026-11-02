#include "cp-algo/math/poly/transform.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <ranges>
#include <span>
#include <sstream>
#include <utility>
#include <vector>
using namespace cp_algo::math;
using T = modint<998244353>;
using P = poly_t<T>;

static_assert(std::ranges::view<decltype(dilate<2>(std::vector<T>{}))>);
using V = decltype(dilate<2>(std::vector<T>{}));
static_assert(std::ranges::random_access_range<V const> && std::ranges::sized_range<V const>);
static_assert(std::ranges::common_range<V> && !std::ranges::contiguous_range<V>);
static_assert(!std::ranges::borrowed_range<V> && !std::ranges::output_range<V, T>);

P oracle(P p, auto&& g) {
    P result;
    for(T c: g) {result += p * c; p = deriv(std::move(p));}
    return result;
}
int main() {
    for(size_t n = 0; n < 12; ++n) for(size_t step: {1, 2, 3, 7, 64}) {
        std::vector<T> a(n);
        for(size_t i = 0; i < n; ++i) {a[i] = i + 1;}
        auto view = dilate(a, step);
        std::vector<T> expected(n ? (n - 1) * step + 1 : 0);
        for(size_t i = 0; i < n; ++i) {expected[i * step] = a[i];}
        assert(view.size() == expected.size() && std::ranges::equal(view, expected));
        assert(view.end() - view.begin() == ptrdiff_t(expected.size()));
        for(size_t i = 0; i < expected.size(); ++i) {
            assert(view[i] == expected[i] && view.begin()[i] == expected[i]);
            assert(*(view.end() - ptrdiff_t(expected.size() - i)) == expected[i]);
        }
        assert(std::ranges::equal(view | std::views::reverse, expected | std::views::reverse));
        if(n) {a[0] = 17; assert(view.front() == T(17)); a[0] = 1;}
        auto owned = dilate(std::vector<T>(a), step);
        auto moved = std::move(owned);
        assert(std::ranges::equal(moved, expected));
        auto [values, offset, stride] = coefficient_layout(std::move(moved));
        assert(stride == step && std::ranges::equal(values, a));
        auto nested = dilate<2>(dilate<3>(std::vector<T>(a)));
        assert(std::ranges::equal(nested, dilate(a, 6)));
    }
    P p(P::Vector{2, 3, 5, 7, 11, 13, 17, 19});
    std::array<T, 3> a{1, 2, 3};
    auto check = [&](auto&& g) {
        auto expected = oracle(p, g);
        P actual = apply_diff(p, g);
        assert(actual == expected);
        auto integrated = integr(apply_diff(p, g));
        assert(integrated[0] == T(0) && deriv(integrated) == expected);
    };
    check(a);
    check(std::span(a));
    struct bare_view: std::ranges::view_base {
        std::span<T> data;
        explicit bare_view(std::span<T> data): data(data) {}
        auto begin() const {return data.begin();}
        auto end() const {return data.end();}
        auto size() const {return data.size();}
    };
    check(dilate<2>(bare_view{std::span(a)}));
    check(std::span<T>{});
    check(dilate<2>(a));
    auto d = dilate<3>(a);
    check(d | std::views::take(5));
    check(d | std::views::drop(1));
    check(d | std::views::reverse);
    check(d | std::views::transform([](T x) {return x + 1;}));
    check(a | std::views::take_while([](T x) {return x != T(3);}));
    auto [fallback, offset, stride] = coefficient_layout(d | std::views::drop(1));
    static_assert(decltype(stride)::value == 1);
    assert(std::ranges::equal(fallback, d | std::views::drop(1)));
    auto [owned, origin, step] = coefficient_layout(dilate<3>(std::vector<T>{4, 5}));
    static_assert(decltype(step)::value == 3);
    assert(owned[0] == T(4) && owned[1] == T(5));
    auto delayed = apply_diff(p, dilate<2>(std::vector<T>{1, 2, 3}));
    assert(P(delayed) == oracle(p, dilate<2>(a)));
    std::istringstream input("1 2 3");
    P from_stream = apply_diff(p, std::ranges::istream_view<int>(input) | std::views::transform([](int x) {return T(x);}));
    assert(from_stream == oracle(p, a));
}
