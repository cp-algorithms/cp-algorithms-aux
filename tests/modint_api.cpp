#include "cp-algo/number_theory/modint.hpp"
#include <cassert>
#include <concepts>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>

using namespace cp_algo::math;

template<class T>
concept fixed_modulus = requires { typename std::integral_constant<typename T::Int, T::mod()>; };
template<class T>
concept switchable = requires { T::switch_mod(7); T::with_mod(7, [] {return 0;}); };
template<class T>
concept montgomery = requires { T::imod(); T::pw128(); T::m_reduce(1); T::m_transform(1); };

static_assert(std::same_as<dynamic_modint<>, modint<0>>);
static_assert(std::same_as<dynamic_modint<long long>, modint<0LL>>);
static_assert(std::same_as<dynamic_modint<int64_t>, modint<int64_t(0)>>);
static_assert(fixed_modulus<modint<17>> && !fixed_modulus<modint<0>>);
static_assert(switchable<modint<0>> && !switchable<modint<17>>);
static_assert(montgomery<modint<0>> && !montgomery<modint<17>>);
static_assert(modint_type<modint<17>> && modint_type<modint<0>> && !modint_type<int>);
static_assert(sizeof(modint<17>) == sizeof(int) && sizeof(modint<0>) == sizeof(int));
static_assert(sizeof(modint<17LL>) == sizeof(long long) && sizeof(modint<0LL>) == sizeof(long long));
static_assert(std::is_trivially_copyable_v<modint<17>> && std::is_trivially_copyable_v<modint<0>>);
static_assert(std::is_standard_layout_v<modint<17>> && std::is_standard_layout_v<modint<0>>);

template<class T>
void arithmetic() {
    T zero{};
    assert(zero.getr() == 0);
    for(int a = 0; a < 17; a++) {
        for(int b = 0; b < 17; b++) {
            T x = a, y = b;
            assert((x + y).getr() == (a + b) % 17);
            assert((x - y).getr() == (a - b + 17) % 17);
            assert((x * y).getr() == a * b % 17);
            assert((-x).getr() == (17 - a) % 17);
            if(b) {assert(x / y * y == x);}
        }
    }
}

int main() {
    arithmetic<modint<17>>();
    arithmetic<modint<17LL>>();
    using D = dynamic_modint<>;
    D::switch_mod(17);
    arithmetic<modint<0>>();
    dynamic_modint<long long>::with_mod(17, [] {arithmetic<modint<0LL>>();});

    assert(D::with_mod(19, [] {return (D(7) * D(7)).getr();}) == 11);
    assert(D::mod() == 17);
    try {
        D::with_mod(23, [] {D::with_mod(29, [] {throw std::runtime_error("test");});});
    } catch(std::runtime_error const&) {}
    assert(D::mod() == 17);

    std::istringstream input("15");
    D value;
    input >> value;
    std::ostringstream output;
    output << value;
    assert(output.str() == "15");

    dynamic_modint<long long>::switch_mod(31);
    assert(D::mod() == 17 && dynamic_modint<long long>::mod() == 31);
    std::thread thread([] {
        assert(D::mod() == 1);
        D::switch_mod(97);
        assert(D(100).getr() == 3);
    });
    thread.join();
    assert(D::mod() == 17);
}
