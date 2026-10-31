#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/util/io.hpp"
#include <bit>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

using namespace cp_algo::math;
using cp_algo::io::read;
using cp_algo::io::read_reduced;

template<class T> uint64_t reference(int64_t value) {
    __int128 r = __int128(value) % __int128(T::mod());
    return uint64_t(r < 0 ? r + T::mod() : r);
}

template<class T, class Int> void explicit_width() {
    using UInt = std::make_unsigned_t<Int>;
    std::vector<Int> values{0, 1, std::numeric_limits<Int>::min(), std::numeric_limits<Int>::max()};
    if constexpr(std::is_signed_v<Int>) {values.push_back(-1);}
    std::mt19937_64 random(47192);
    for(int i = 0; i < 128; i++) {values.push_back(std::bit_cast<Int>(UInt(random())));}
    std::stringstream input;
    for(Int value: values) {input << value << ' ';}
    for(Int value: values) {
        T x;
        read<Int>(x, input);
        __int128 expected = __int128(value) % __int128(T::mod());
        if(expected < 0) {expected += T::mod();}
        assert(input && x.getr() == uint64_t(expected));
        assert((x * x).getr() == __uint128_t(expected) * expected % uint64_t(T::mod()));
    }
}

template<class T> void check() {
    explicit_width<T, int32_t>();
    explicit_width<T, uint32_t>();
    explicit_width<T, int64_t>();
    explicit_width<T, uint64_t>();
    using Input = std::make_signed_t<typename T::UInt>;
    constexpr int64_t lo = std::numeric_limits<Input>::min();
    constexpr int64_t hi = std::numeric_limits<Input>::max();
    int64_t m = T::mod();
    {
        std::stringstream input;
        for(int64_t v : {int64_t(0), m / 2, m - 1}) {input << v << ' ';}
        for(int64_t v : {int64_t(0), m / 2, m - 1}) {
            T value;
            read_reduced(value, input);
            assert(input && value.getr() == uint64_t(v));
        }
        for(std::string text : {"", "word", "18446744073709551616"}) {
            std::istringstream bad(text);
            T value = 7;
            auto before = value.getr();
            read_reduced(value, bad);
            assert(bad.fail() && value.getr() == before);
        }
    }
    std::vector<int64_t> values{0, 1, -1, m - 1, m, -m, -m + 1, lo, hi};
    for(int64_t v : {int64_t(INT32_MIN), int64_t(INT32_MAX), int64_t(3000000000), int64_t(UINT32_MAX),
                    int64_t(UINT32_MAX) + 1, -int64_t(UINT32_MAX) - 1}) {
        if(lo <= v && v <= hi) {values.push_back(v);}
    }
    if(m < hi) {values.push_back(m + 1);}
    std::mt19937_64 random(627491);
    for(int i = 0; i < 1024; i++) {
        values.push_back(std::bit_cast<Input>(typename T::UInt(random())));
    }
    std::stringstream input;
    for(auto v : values) { input << v << ' '; }
    for(auto v : values) {
        T value = 0;
        read<Input>(value, input);
        assert(input);
        auto expected = reference<T>(v);
        assert(value.getr() == expected);
        assert((value + T(1)).getr() == (expected + 1) % uint64_t(m));
        assert((value * value).getr() == __uint128_t(expected) * expected % uint64_t(m));
    }

    std::vector<std::string> invalid{"", " \t\n", "word", "+", "-", "9223372036854775808", "-9223372036854775809"};
    if constexpr(sizeof(Input) == 4) {
        invalid.insert(invalid.end(), {"2147483648", "-2147483649", "3000000000"});
    }
    for(auto const& text : invalid) {
        std::istringstream bad(text);
        T value = 7;
        auto before = value.getr();
        read<Input>(value, bad);
        assert(bad.fail());
        assert(value.getr() == before);
    }
    for(auto state : {std::ios::failbit, std::ios::badbit, std::ios::eofbit}) {
        std::istringstream bad("123");
        bad.setstate(state);
        T value = 7;
        auto before = value.getr();
        read<Input>(value, bad);
        assert(bad.fail() && value.getr() == before);
    }
    {
        std::istringstream bad("invalid");
        bad.exceptions(std::ios::failbit);
        T value = 7;
        auto before = value.getr();
        bool threw = false;
        try { read<Input>(value, bad); } catch(std::ios_base::failure const&) { threw = true; }
        assert(threw && value.getr() == before);
    }
    {
        std::istringstream input("  +15 -0 12x");
        T a, b, c;
        read<Input>(a, input); read<Input>(b, input); read<Input>(c, input);
        assert(a.getr() == reference<T>(15) && b.getr() == 0 && c.getr() == reference<T>(12));
        assert(input.peek() == 'x');
        auto before = c.getr();
        read<Input>(c, input);
        assert(input.fail() && c.getr() == before);
    }
    {
        std::istringstream input(sizeof(Input) == 4 ? "ff -ff 7fffffff -80000000" :
                                                    "ff -ff 7fffffffffffffff -8000000000000000");
        input >> std::hex;
        for(int64_t v : {int64_t(255), int64_t(-255), hi, lo}) {
            T value;
            read<Input>(value, input);
            assert(!input.fail() && value.getr() == reference<T>(v));
        }
        assert(input.eof());
    }
}

int main() {
    check<modint<1>>();
    check<modint<2>>();
    check<modint<998244353>>();
    check<modint<998244353u>>();
    check<modint<INT32_MAX>>();
    check<modint<2305843009213693951LL>>();
    check<modint<2305843009213693951ULL>>();
    check<modint<INT64_MAX>>();
    for(int m : {1, 2, 3, 10, 998244353, 1073741823, 1073741825, 1879048201, INT32_MAX}) {
        dynamic_modint<>::with_mod(m, [] { check<dynamic_modint<>>(); });
    }
    for(int64_t m : {int64_t(1), int64_t(2), int64_t(998244353), int64_t(2305843009213693951), INT64_MAX - 1, INT64_MAX}) {
        dynamic_modint<int64_t>::with_mod(m, [] { check<dynamic_modint<int64_t>>(); });
    }
}
