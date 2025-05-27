#ifndef CP_ALGO_STRUCTURES_BITPACK_HPP
#define CP_ALGO_STRUCTURES_BITPACK_HPP
#include "../structures/bit_array.hpp"
#include <cstdint>
#include <cstddef>
#include <string>
#include <array>
namespace cp_algo::structures {
    template<size_t n, typename Int = uint64_t>
    struct bitpack: bit_array<n, Int> {
        using Base = bit_array<n, Int>;
        using Base::width, Base::blocks, Base::data;
        auto operator <=> (bitpack const& t) const = default;

        bitpack() {}
        bitpack(std::string bits) {
            size_t rem = size(bits) % width;
            if(rem) {
                bits += std::string(width - rem, '0');
            }
            for(size_t i = 0, pos = 0; pos < size(bits); i++, pos += width) {
                for(size_t j = width; j; j--) {
                    data[i] *= 2;
                    data[i] ^= bits[pos + j - 1] == '1';
                }
            }
        }

        bitpack& xor_hint(bitpack const& t, size_t hint) {
            for(size_t i = hint / width; i < blocks; i++) {
                data[i] ^= t.data[i];
            }
            return *this;
        }
        bitpack& operator ^= (bitpack const& t) {
            return xor_hint(t, 0);
        }
        bitpack operator ^ (bitpack const& t) const {
            return bitpack(*this) ^= t;
        }

        std::string to_string() const {
            std::string res(blocks * width, '0');
            for(size_t i = 0, pos = 0; i < blocks; i++, pos += width) {
                Int block = data[i];
                for(size_t j = 0; j < width; j++) {
                    res[pos + j] = '0' + block % 2;
                    block /= 2;
                }
            }
            res.resize(n);
            return res;
        }

        size_t ctz() const {
            size_t res = 0;
            size_t i = 0;
            while(i < blocks && data[i] == 0) {
                res += width;
                i++;
            }
            if(i < blocks) {
                res += std::countr_zero(data[i]);
            }
            return std::min(res, n);
        }
    };
}
#endif // CP_ALGO_STRUCTURES_BITPACK_HPP
