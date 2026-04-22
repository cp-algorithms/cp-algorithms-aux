#ifndef CP_ALGO_UTIL_GENERATOR_HPP
#define CP_ALGO_UTIL_GENERATOR_HPP

#include "big_alloc.hpp"
#include <generator>

namespace cp_algo {
    template<typename Ref, typename V = void>
    using big_generator = std::generator<Ref, V, big_alloc<std::byte>>;
}

namespace std::ranges {
    template<typename Ref, typename V>
    elements_of(cp_algo::big_generator<Ref, V>&&) -> elements_of<cp_algo::big_generator<Ref, V>&&, cp_algo::big_alloc<std::byte>>;
}

#endif // CP_ALGO_UTIL_GENERATOR_HPP
