#ifndef CP_ALGO_UTIL_IO_HPP
#define CP_ALGO_UTIL_IO_HPP
#include <iostream>
#include <cstddef>
#include <ranges>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

namespace cp_algo::io {
    namespace impl {
        // The extra type tag keeps ADL customizations distinct from the public
        // entry points, even when those are imported into a stream's namespace.
        void read() = delete;
        void read_reduced() = delete;
        void write() = delete;

        template<class Int, bool reduced, class T, class Stream>
        constexpr bool custom_reader = [] {
            if constexpr(reduced) {
                return requires(T& x, Stream& in) {read_reduced(x, in, std::type_identity<void>{});};
            } else {
                return requires(T& x, Stream& in) {read(x, in, std::type_identity<Int>{});};
            }
        }();

        template<class Int, bool reduced, class T, class Stream>
        consteval bool readable() {
            if constexpr(custom_reader<Int, reduced, T, Stream>) {
                return true;
            } else if constexpr(requires(T& x, std::istream& in) {in >> x;}) {
                if constexpr(!std::is_void_v<Int> && std::is_arithmetic_v<T>) {
                    return std::is_arithmetic_v<Int> && requires(Int& value, T& x, std::istream& in) {
                        in >> value;
                        x = static_cast<T>(value);
                    };
                } else {return true;}
            } else if constexpr(std::ranges::input_range<T>) {
                return readable<Int, reduced, std::remove_reference_t<std::ranges::range_reference_t<T>>, Stream>();
            } else if constexpr(requires {std::tuple_size<std::remove_cv_t<T>>::value;}) {
                return []<size_t... I>(std::index_sequence<I...>) {
                    return (readable<Int, reduced, std::remove_reference_t<decltype(std::get<I>(std::declval<T&>()))>, Stream>() && ...);
                }(std::make_index_sequence<std::tuple_size_v<std::remove_cv_t<T>>>{});
            } else {return false;}
        }

        template<class Int, bool reduced, class T, class Stream>
        void read_value(T& x, Stream& in) {
            if constexpr(custom_reader<Int, reduced, T, Stream>) {
                if constexpr(reduced) {read_reduced(x, in, std::type_identity<void>{});}
                else {read(x, in, std::type_identity<Int>{});}
            } else if constexpr(requires(std::istream& native) {native >> x;}) {
                if constexpr(!std::is_void_v<Int> && std::is_arithmetic_v<T>) {
                    Int value;
                    if(in >> value) {x = static_cast<T>(value);}
                } else {in >> x;}
            } else if constexpr(std::ranges::input_range<T>) {
                for(auto& value: x) {read_value<Int, reduced>(value, in);}
            } else {
                std::apply([&](auto&... value) {(read_value<Int, reduced>(value, in), ...);}, x);
            }
        }

        template<class T, class Stream>
        constexpr bool custom_writer = requires(T const& x, Stream& out) {write(x, out, std::type_identity<void>{});};

        // Numeric C arrays also decay to printable pointers; visit their elements.
        template<class T>
        constexpr bool native_writer = requires(T const& x, std::ostream& out) {out << x;} &&
            (!std::is_array_v<T> || std::is_convertible_v<T const&, std::string_view>);

        template<class T, class Stream>
        consteval bool writable() {
            if constexpr(custom_writer<T, Stream> || native_writer<T>) {
                return true;
            } else if constexpr(std::ranges::input_range<T const>) {
                return writable<std::remove_reference_t<std::ranges::range_reference_t<T const>>, Stream>();
            } else if constexpr(requires {std::tuple_size<T>::value;}) {
                return []<size_t... I>(std::index_sequence<I...>) {
                    return (writable<std::tuple_element_t<I, T>, Stream>() && ...);
                }(std::make_index_sequence<std::tuple_size_v<T>>{});
            } else {return false;}
        }

        template<class T, class Stream>
        void write_value(T const& x, Stream& out, bool& first) {
            if constexpr(custom_writer<T, Stream> || native_writer<T>) {
                if(!first) {out << ' ';}
                if constexpr(custom_writer<T, Stream>) {write(x, out, std::type_identity<void>{});}
                else {out << x;}
                first = false;
            } else if constexpr(std::ranges::input_range<T const>) {
                for(auto const& value: x) {write_value(value, out, first);}
            } else {
                std::apply([&](auto const&... value) {(write_value(value, out, first), ...);}, x);
            }
        }
    }

    // Containers must be sized already. Int selects numeric input width; strings
    // keep their ordinary stream extraction. Ordinary numeric results must fit
    // their destination type. With no Int, use each leaf's reader.
    template<class Int = void, class T, class Stream = decltype(std::cin)>
        requires(impl::readable<Int, false, T, Stream>())
    void read(T& x, Stream& in = std::cin) {
        impl::read_value<Int, false>(x, in);
    }

    // Modular leaves must already be in [0, mod()); ordinary leaves use their reader.
    template<class T, class Stream = decltype(std::cin)>
        requires(impl::readable<void, true, T, Stream>())
    void read_reduced(T& x, Stream& in = std::cin) {
        impl::read_value<void, true>(x, in);
    }

    // Flatten pairs, tuples and ranges, separating leaves by spaces. No newline.
    template<class T, class Stream = decltype(std::cout)>
        requires(impl::writable<T, Stream>())
    void write(T const& x, Stream& out = std::cout) {
        bool first = true;
        impl::write_value(x, out, first);
    }
}
#endif // CP_ALGO_UTIL_IO_HPP
