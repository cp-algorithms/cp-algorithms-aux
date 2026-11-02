#ifndef CP_ALGO_MATH_POLY_COEFFICIENT_VIEW_HPP
#define CP_ALGO_MATH_POLY_COEFFICIENT_VIEW_HPP
#include <cassert>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <ranges>
#include <type_traits>
#include <utility>

namespace cp_algo::math {
namespace poly::impl {
    template<std::ranges::range R> requires(std::is_lvalue_reference_v<R> || std::ranges::viewable_range<R>)
    auto range_view(R&& r) {
        if constexpr(std::ranges::viewable_range<R>) {return std::views::all(std::forward<R>(r));}
        else {return std::ranges::ref_view(r);}
    }
}
    // Possible support at offset + j*stride; stored values may also be zero.
    template<std::ranges::view V, typename Stride>
    struct coefficient_layout_t {
        V values;
        size_t offset;
        Stride stride;
        size_t extent(size_t n) const {
            if(!n) {return 0;}
            size_t limit = size_t(std::numeric_limits<ptrdiff_t>::max());
            assert(stride > 0 && offset < limit);
            assert(n - 1 <= (limit - 1 - offset) / stride);
            return offset + (n - 1) * stride + 1;
        }
    };
    template<std::ranges::range R> requires(std::is_lvalue_reference_v<R> || std::ranges::viewable_range<R>)
    auto coefficient_layout(R&& r) {
        if constexpr(requires { std::forward<R>(r).coefficient_layout(); }) {
            return std::forward<R>(r).coefficient_layout();
        } else {
            return coefficient_layout_t{poly::impl::range_view(std::forward<R>(r)), 0,
                                        std::integral_constant<size_t, 1>{}};
        }
    }

    // A read-only view of x^offset * p(x^stride), including the zeros.
    // Store ranges as std::views::all does; also borrow move-only lvalue views.
    // Iterators refer to this view and must not outlive it or its underlying storage.
    template<std::ranges::view V, typename Stride>
        requires std::ranges::random_access_range<V const> && std::ranges::sized_range<V const>
    struct coefficient_view: std::ranges::view_interface<coefficient_view<V, Stride>> {
        V coefficients;
        size_t offset;
        Stride stride;
        using value_type = std::ranges::range_value_t<V>;

        coefficient_view(V coefficients, size_t offset, Stride stride):
            coefficients(std::move(coefficients)), offset(offset), stride(stride) {
            assert(stride > 0);
        }
        size_t size() const {return coefficient_layout().extent(std::ranges::size(coefficients));}
        value_type operator[](size_t i) const {
            if(i < offset || (i - offset) % stride) {return {};}
            return std::ranges::begin(coefficients)[(i - offset) / stride];
        }

        struct iterator {
            using value_type = coefficient_view::value_type;
            using difference_type = ptrdiff_t;
            using iterator_concept = std::random_access_iterator_tag;
            using iterator_category = std::input_iterator_tag;
            coefficient_view const* view = nullptr;
            difference_type index = 0;
            value_type operator*() const {return (*view)[index];}
            value_type operator[](difference_type n) const {return (*view)[index + n];}
            iterator& operator+=(difference_type n) {index += n; return *this;}
            iterator& operator-=(difference_type n) {index -= n; return *this;}
            iterator& operator++() {return *this += 1;}
            iterator& operator--() {return *this -= 1;}
            iterator operator++(int) {auto old = *this; ++*this; return old;}
            iterator operator--(int) {auto old = *this; --*this; return old;}
            friend iterator operator+(iterator p, difference_type n) {return p += n;}
            friend iterator operator+(difference_type n, iterator p) {return p += n;}
            friend iterator operator-(iterator p, difference_type n) {return p -= n;}
            friend difference_type operator-(iterator a, iterator b) {return a.index - b.index;}
            auto operator<=>(iterator const&) const = default;
        };
        iterator begin() const {return {this, 0};}
        iterator end() const {return {this, ptrdiff_t(size())};}

        auto coefficient_layout(this auto&& self) {
            return coefficient_layout_t{poly::impl::range_view(std::forward<decltype(self)>(self).coefficients), self.offset, self.stride};
        }
    };

namespace poly::impl {
    template<typename P>
    auto layout_of(P&& p) {
        if constexpr(std::ranges::range<P>) {return coefficient_layout(std::forward<P>(p));}
        else {return coefficient_layout(std::forward<P>(p).a);}
    }
    template<size_t N, typename X>
    auto scale_stride(X x) {
        assert(size_t(x) <= SIZE_MAX / N);
        if constexpr(requires { X::value; }) {return std::integral_constant<size_t, N * X::value>{};}
        else {return N * size_t(x);}
    }
}
    template<typename P>
    auto dilate(P&& p, size_t stride) {
        auto layout = poly::impl::layout_of(std::forward<P>(p));
        assert(stride > 0 && layout.offset <= SIZE_MAX / stride && size_t(layout.stride) <= SIZE_MAX / stride);
        return coefficient_view{std::move(layout.values), layout.offset * stride, size_t(layout.stride) * stride};
    }
    template<size_t Stride, typename P>
    auto dilate(P&& p) {
        static_assert(Stride > 0);
        auto layout = poly::impl::layout_of(std::forward<P>(p));
        assert(layout.offset <= SIZE_MAX / Stride);
        return coefficient_view{std::move(layout.values), layout.offset * Stride, poly::impl::scale_stride<Stride>(layout.stride)};
    }
    // View of x^k * p(x); an empty coefficient range remains empty.
    template<typename P>
    auto mul_xk(P&& p, size_t k) {
        auto layout = poly::impl::layout_of(std::forward<P>(p));
        assert(k <= SIZE_MAX - layout.offset);
        return coefficient_view{std::move(layout.values), layout.offset + k, layout.stride};
    }
}
#endif // CP_ALGO_MATH_POLY_COEFFICIENT_VIEW_HPP
