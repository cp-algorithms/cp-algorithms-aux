#include "cp-algo/util/io.hpp"
#include "cp-algo/number_theory/modint.hpp"
#include <array>
#include <cassert>
#include <cstdint>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using mint = cp_algo::math::modint<17>;
using cp_algo::io::read;
using cp_algo::io::read_reduced;
using cp_algo::io::write;

template<class T> concept implicit_reader = requires(T& x) {cp_algo::io::read(x);};
template<class T> concept stream_reader = requires(T& x, std::istream& in) {in >> x;};
static_assert(implicit_reader<int> && implicit_reader<std::vector<int>>);
static_assert(!implicit_reader<mint> && !implicit_reader<std::vector<mint>>);
static_assert(!implicit_reader<std::pair<int, mint>> && !stream_reader<mint>);
static_assert(!implicit_reader<std::map<int, int>>);

// An unconstrained stream operator must not make a container look like a scalar.
struct input_backend {
    std::istringstream input;
    template<class T> input_backend& operator>>(T& x) {input >> x; return *this;}
    explicit operator bool() const {return bool(input);}
};
struct output_backend {
    std::ostringstream output;
    template<class T> output_backend& operator<<(T const& x) {output << x; return *this;}
};

template<class In, class Out>
void nested(In& input, Out& output) {
    std::tuple<int, std::string, std::vector<mint>> values{0, "", std::vector<mint>(3)};
    read<int64_t>(values, input);
    assert(std::get<0>(values) == 7 && std::get<1>(values) == "hello");
    assert(std::get<2>(values)[0].getr() == 8);
    assert(std::get<2>(values)[1].getr() == 8);
    assert(std::get<2>(values)[2].getr() == 9);
    write(values, output);
}

int main() {
    std::string text = "7 hello -9223372036854775808 9223372036854775807 2147483648";
    std::istringstream input(text);
    std::ostringstream output;
    nested(input, output);
    assert(output.str() == "7 hello 8 8 9");
    input_backend backend{std::istringstream(text)};
    output_backend sink;
    nested(backend, sink);
    assert(sink.output.str() == output.str());

    std::array<std::pair<int, mint>, 2> pairs;
    std::istringstream reduced("3 15 4 16");
    read_reduced(pairs, reduced);
    assert(pairs[0].first == 3 && pairs[0].second.getr() == 15);
    assert(pairs[1].first == 4 && pairs[1].second.getr() == 16);

    std::vector<mint> wide(2);
    std::istringstream unsigned_input("18446744073709551615 4294967295");
    read<uint64_t>(wide, unsigned_input);
    assert(wide[0].getr() == 0 && wide[1].getr() == 0);

    int ordinary[2][2]{};
    std::istringstream ordinary_input("1 2 3 4");
    read(ordinary, ordinary_input);
    std::ostringstream ordinary_output;
    write(ordinary, ordinary_output);
    assert(ordinary_output.str() == "1 2 3 4");

    std::vector<mint> partial{5, 7};
    std::istringstream bad("2 invalid");
    read<int>(partial, bad);
    assert(bad.fail() && partial[0].getr() == 2 && partial[1].getr() == 7);

    std::ostringstream empty;
    write(std::tuple{}, empty);
    write(std::vector<int>{}, empty);
    assert(empty.str().empty());

    std::istringstream defaults("-1 16");
    auto* previous = std::cin.rdbuf(defaults.rdbuf());
    mint scalar;
    read<int>(scalar);
    assert(scalar.getr() == 16);
    read_reduced(scalar);
    assert(scalar.getr() == 16);
    std::cin.rdbuf(previous);
    std::cin.clear();
}
