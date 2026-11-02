#include "cp-algo/math/fps.hpp"
#include "cp-algo/math/laurent.hpp"
#include "cp-algo/math/poly/transform.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <utility>
using namespace cp_algo::math;
using T = modint<998244353>;
using P = poly_t<T>;
using V = decltype(std::declval<fps<T>>().prefix_view(0));
static_assert(std::ranges::view<V> && std::ranges::random_access_range<V const>);
static_assert(std::ranges::sized_range<V> && !std::ranges::borrowed_range<V>);
static_assert(!std::ranges::output_range<V, T>);

P oracle(P p, auto const& g) {
    P answer;
    for(T x: g) {answer += p * x; p = deriv(std::move(p));}
    return answer;
}

int main() {
    size_t calls = 0, allowed = 0;
    fps<T> f([&](size_t i, auto const&) {
        assert(i == calls++ && i <= allowed);
        return T(i + 1);
    });
    auto empty = f.prefix_view(0);
    auto view = f.prefix_view(100);
    assert(empty.empty() && view.size() == 100 && calls == 0);
    assert(view[0] == T(1) && calls == 1);
    allowed = 3;
    assert(view[3] == T(4) && calls == 4);
    auto copy = view;
    f = fps<T>();
    allowed = 1000;
    assert(copy[99] == T(100) && view[99] == T(100) && calls == 100);
    assert(std::ranges::equal(view | std::views::take(4), std::array<T, 4>{1, 2, 3, 4}));
    assert(view.end() - view.begin() == 100);

    size_t generated = 0;
    fps<T> sparse([&](size_t i, auto const&) {++generated; return T(i + 1);});
    auto spaced = mul_xk(dilate<2>(sparse.prefix_view(3)), 1000000);
    assert(spaced.size() == 1000005 && generated == 0);
    assert(spaced[999999] == T(0) && spaced[1000001] == T(0) && generated == 0);
    assert(spaced[1000004] == T(3) && generated == 3);
    assert(std::ranges::size(coefficient_layout(spaced).values) == 3);

    P aliased(P::Vector{2, 3, 5, 7});
    P wanted = aliased + aliased.mul_xk(3);
    fps<T> source([&](size_t i, auto const&) {return aliased[int(i)];});
    aliased += mul_xk(source.prefix_view(4), 3);
    assert(aliased == wanted);

    auto temporary = fps<int>([p = std::make_unique<int>(7)](size_t i, auto const&) {return *p + int(i);}).prefix_view(3);
    assert(std::ranges::equal(temporary, std::array{7, 8, 9}));
    auto growing = fps<T>([](size_t i, auto const&) {return T(i + 2);});
    auto prefix = growing.prefix_view(4);
    assert(prefix[0] == T(2));
    assert(growing[100000] == T(100002));
    assert(std::ranges::equal(prefix, std::array<T, 4>{2, 3, 4, 5}));

    auto fixed = fps<T>(P(P::Vector{2, 3})).prefix_view(5);
    assert(std::ranges::equal(fixed, std::array<T, 5>{2, 3, 0, 0, 0}));
    auto layout = coefficient_layout(fixed);
    assert(layout.offset == 0 && layout.stride == 1 && std::ranges::size(layout.values) == 5);

    P p(P::Vector{2, 3, 5, 7, 11, 13, 17, 19});
    for(size_t n = 0; n <= 10; ++n) for(size_t stride: {1, 2, 3}) for(size_t offset: {0, 1, 3, 12}) {
        auto series = inv(fps<T>(P(P::Vector{1, -1})));
        auto g = mul_xk(dilate(series.prefix_view(n), stride), offset);
        P expected = oracle(p, g);
        assert(P(apply_diff(p, g)) == expected);
        assert(deriv(integr(apply_diff(p, g))) == expected);
        auto q = p;
        q += g;
        assert(q == p + P(std::ranges::to<P::Vector>(g)));
    }
    auto h = mul_xk(dilate<2>(growing.prefix_view(3)), 1);
    static_assert(decltype(coefficient_layout(h).stride)::value == 2);
    assert(std::ranges::equal(h, std::array<T, 6>{0, 2, 0, 3, 0, 4}));
    laurent<T> shifted{growing, 3};
    auto aligned = mul_xk(shifted.series.prefix_view(4), size_t(shifted.offset));
    for(size_t i = 0; i < aligned.size(); ++i) {assert(aligned[i] == shifted[int64_t(i)]);}

    fps<T> fibonacci([&](size_t i, auto const& previous) {
        return i < 2 ? T(1) : previous[i - 1] + previous[i - 2];
    });
    assert(std::ranges::equal(fibonacci.prefix_view(6), std::array<T, 6>{1, 1, 2, 3, 5, 8}));
    fps<T> cycle;
    cycle = fps<T>([&](size_t i, auto const&) {return cycle.prefix_view(i + 2)[i + 1];});
    for(int attempt = 0; attempt < 2; ++attempt) {
        bool rejected = false;
        try {(void)cycle.prefix_view(1)[0];} catch(std::logic_error const&) {rejected = true;}
        assert(rejected);
    }
    size_t attempts = 0;
    fps<T> retry([&](size_t i, auto const&) {
        if(i == 3 && attempts++ == 0) {throw std::runtime_error("retry");}
        return T(i);
    });
    auto retried = retry.prefix_view(6);
    bool rejected = false;
    try {(void)retried[5];} catch(std::runtime_error const&) {rejected = true;}
    assert(rejected && retried[2] == T(2) && retried[5] == T(5));
}
