#include "cp-algo/math/poly/coefficient_view.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <ranges>
#include <span>
#include <utility>
#include <vector>
using namespace cp_algo::math;

template<typename R>
concept has_layout = requires(R&& r) {coefficient_layout(std::forward<R>(r));};
using Owned = decltype(mul_xk(dilate<2>(std::vector<int>{}), 1));
using Borrowed = decltype(dilate<2>(std::span<int>{}));
static_assert(std::ranges::view<Owned> && std::movable<Owned> && !std::copyable<Owned>);
static_assert(std::ranges::random_access_range<Owned const> && std::ranges::sized_range<Owned const>);
static_assert(std::ranges::common_range<Owned> && !std::ranges::contiguous_range<Owned>);
static_assert(!std::ranges::borrowed_range<Owned> && !std::ranges::borrowed_range<Borrowed>);
static_assert(!std::ranges::output_range<Borrowed, int>);
static_assert(has_layout<Owned&> && has_layout<Owned const&> && has_layout<Owned>);
static_assert(!has_layout<Owned const> && !has_layout<std::initializer_list<int>>);
static_assert(has_layout<std::initializer_list<int>&>);

auto from_span(std::span<int> s) {return mul_xk(dilate<2>(s), 1);}
auto owned() {return mul_xk(dilate<2>(std::vector<int>{2, 0, 7}), 1);}

int main() {
    std::array<int, 3> a{2, 0, 7};
    auto v = from_span(a);
    auto copy = v;
    auto q = owned();
    auto moved = std::move(q);
    const auto& cv = v;
    assert(std::ranges::equal(cv, std::array{0, 2, 0, 0, 0, 7}));
    assert(std::ranges::equal(v, moved) && std::ranges::equal(copy, moved));
    assert(v.end() - v.begin() == 6 && v.begin() < v.end());
    auto it = v.begin();
    assert(*++it == 2 && *it++ == 2 && *it == 0);
    assert(*--it == 2 && *it-- == 2 && it == v.begin());
    assert(*(3 + it) == 0 && (it + 5)[-4] == 2);
    decltype(it) empty1, empty2;
    assert(empty1 == empty2);
    a[0] = 5;
    assert(v[1] == 5 && copy[1] == 5 && moved[1] == 2);

    auto layout = coefficient_layout(moved);
    static_assert(decltype(layout.stride)::value == 2);
    assert(layout.offset == 1 && std::ranges::size(layout.values) == 3);
    auto nested = dilate<3>(mul_xk(moved, 2));
    auto nested_layout = coefficient_layout(nested);
    static_assert(decltype(nested_layout.stride)::value == 6);
    assert(nested_layout.offset == 9 && nested.size() == 22 && nested[21] == 7);
    auto acquired = coefficient_layout(owned());
    assert(std::ranges::equal(acquired.values, std::array{2, 0, 7}));
    const auto const_owned = owned();
    assert(std::ranges::equal(coefficient_layout(const_owned).values, acquired.values));

    // Ordinary adapters must account for virtual zeros when the support hint is lost.
    auto shifted = v | std::views::drop(1);
    auto dense = coefficient_layout(shifted);
    static_assert(decltype(dense.stride)::value == 1);
    assert(dense.offset == 0 && std::ranges::size(dense.values) == 5);
    assert(std::ranges::equal(dilate<2>(shifted), std::array{5, 0, 0, 0, 0, 0, 0, 0, 7}));
    assert(std::ranges::equal(v | std::views::reverse, std::array{7, 0, 0, 0, 5, 0}));
    assert(std::ranges::equal(v | std::views::transform([](int x) {return x + 1;}),
                              std::array{1, 6, 1, 1, 1, 8}));

    size_t limit = size_t(std::numeric_limits<ptrdiff_t>::max());
    auto high = mul_xk(std::views::single(9), limit - 1);
    assert(high.size() == limit && high.end() - high.begin() == ptrdiff_t(limit));
    assert(high[limit - 1] == 9 && high[limit - 2] == 0);
    auto empty = mul_xk(std::views::empty<int>, SIZE_MAX);
    assert(empty.empty() && empty.begin() == empty.end());
}
