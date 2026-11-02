#include "cp-algo/math/poly/transform.hpp"
#include <array>
#include <cassert>
#include <initializer_list>
#include <iterator>
#include <ranges>
#include <span>
#include <sstream>
#include <utility>
#include <vector>
using namespace cp_algo::math;
using T = modint<998244353>;
using P = poly_t<T>;
template<typename R>
concept has_layout = requires(R&& r) {coefficient_layout(std::forward<R>(r));};
using Owned = std::ranges::owning_view<std::vector<T>>;
static_assert(has_layout<Owned&> && has_layout<Owned const&> && has_layout<Owned>);
static_assert(!has_layout<Owned const>);

P materialize(auto const& v) {return std::ranges::to<P::Vector>(v);}
P differential(P p, auto const& coefficients) {
    P q;
    for(T c: coefficients) {q += p * c; p = deriv(std::move(p));}
    return q;
}
int main() {
    P p(P::Vector{2, 3, 5, 7, 11, 13, 17, 19});
    for(size_t n = 0; n < 12; ++n) for(size_t offset: {0, 1, 3, 23}) for(size_t stride: {1, 2, 3, 7}) {
        std::vector<T> a(n);
        for(size_t i = 0; i < n; ++i) {a[i] = i % 3 ? i + 1 : 0;}
        auto view = mul_xk(dilate(a, stride), offset);
        P::Vector expected(n ? offset + (n - 1) * stride + 1 : 0);
        for(size_t i = 0; i < n; ++i) {expected[offset + i * stride] = a[i];}
        assert(view.size() == expected.size() && std::ranges::equal(view, expected));
        auto layout = coefficient_layout(view);
        assert(layout.offset == offset && layout.stride == stride && std::ranges::equal(layout.values, a));
        assert(std::ranges::equal(dilate(mul_xk(a, offset), stride), mul_xk(dilate(a, stride), offset * stride)));
        auto q = p;
        q += view;
        assert(q == p + P(expected));
        q = p;
        q -= view;
        assert(q == p - P(expected));
        q = P{};
        q -= view;
        assert(q == -P(expected));
        assert(P(apply_diff(p, view)) == differential(p, expected));
        auto integrated = integr(apply_diff(p, view));
        assert(integrated[0] == T(0) && deriv(integrated) == differential(p, expected));
        auto slice = view | std::views::drop(1) | std::views::take(12);
        q = p;
        q += slice;
        assert(q == p + materialize(slice));
    }
    for(size_t offset: {0, 1, 3, 20}) for(size_t stride: {1, 2, 3}) {
        auto q = p;
        P expected = p + materialize(mul_xk(dilate(p, stride), offset));
        q += mul_xk(dilate(q, stride), offset);
        assert(q == expected);
        q = p;
        expected = p - materialize(mul_xk(dilate(p, stride), offset));
        q -= mul_xk(dilate(q, stride), offset);
        assert(q == expected);
        q = p;
        auto tail = std::span(q.a).subspan(1);
        expected = p + materialize(mul_xk(tail, offset));
        q += mul_xk(tail, offset);
        assert(q == expected);
        q = p;
        auto reversed = q.a | std::views::reverse;
        expected = p - materialize(mul_xk(reversed, offset));
        q -= mul_xk(reversed, offset);
        assert(q == expected);
    }
    auto owned = mul_xk(dilate<2>(P(P::Vector{4, 0, 5})), 3);
    auto [values, offset, stride] = coefficient_layout(owned);
    static_assert(decltype(stride)::value == 2);
    assert(offset == 3 && std::ranges::size(values) == 3);
    auto nested = mul_xk(owned, 4);
    auto nested_layout = coefficient_layout(nested);
    static_assert(decltype(nested_layout.stride)::value == 2);
    assert(nested_layout.offset == 7 && std::ranges::size(nested_layout.values) == 3);
    auto q = p;
    q += owned;
    assert(q == p + materialize(owned));
    auto transformed = std::ranges::ref_view(owned) | std::views::transform([](T c) {return c + 1;});
    q = p;
    q += transformed;
    assert(q == p + materialize(transformed));
    auto [fallback, start, step] = coefficient_layout(transformed);
    static_assert(decltype(step)::value == 1);
    assert(start == 0);
    auto empty = mul_xk(std::vector<T>{}, size_t(1'000'000));
    assert(empty.empty());
    q = p;
    q += empty;
    assert(q == p);
    assert(P(apply_diff(p, mul_xk(std::vector<T>{1}, size_t(1'000'000)))).is_zero());
    auto counted = std::ranges::subrange(std::counted_iterator(p.a.data(), 3), std::default_sentinel);
    q = p;
    q += counted;
    assert(q == p + p.mod_xk(3));
    std::istringstream input("1 2 3");
    q = p;
    q += std::ranges::istream_view<int>(input);
    assert(q == p + P(P::Vector{1, 2, 3}));
    q = p;
    q += std::initializer_list<T>{1, 2, 3};
    assert(q == p + P(P::Vector{1, 2, 3}));
    const auto fixed = mul_xk(dilate<2>(std::vector<T>{1, 2, 3}), 1);
    auto expected = materialize(fixed);
    q = p;
    q += std::move(fixed);
    assert(q == p + expected);
    assert(P(apply_diff(p, std::move(fixed))) == differential(p, expected.a));
    assert(P(apply_diff(p, std::initializer_list<T>{1, 2, 3})) == differential(p, std::array<T, 3>{1, 2, 3}));
    struct mutable_view: std::ranges::view_base {
        std::span<T> data;
        explicit mutable_view(std::span<T> data): data(data) {}
        auto begin() {return data.begin();}
        auto end() {return data.end();}
        auto size() {return data.size();}
    };
    q = p;
    q += mutable_view(std::span(p.a).first(3));
    assert(q == p + p.mod_xk(3));
    q = p;
    q += mutable_view(std::span(q.a).first(3));
    assert(q == p + p.mod_xk(3));
    assert(P(apply_diff(p, mutable_view(std::span(p.a).first(3)))) == differential(p, p.mod_xk(3).a));
}
