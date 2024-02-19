
#ifndef CP_ALGO_DATA_STRUCTURES_BITPACK_HPP
#define CP_ALGO_DATA_STRUCTURES_BITPACK_HPP
#include <cstdint>
#include <cstddef>
#include <string>
#include <array>
#include <bit>
namespace cp_algo::data_structures {
    template<size_t n, typename Int = uint64_t>
    struct bitpack {
        static constexpr uint8_t bits_per_block = 8 * sizeof(Int);
        static constexpr uint32_t blocks = (n + bits_per_block - 1) / bits_per_block;
        std::array<Int, blocks> data;

        auto operator <=> (bitpack const& t) const = default;

        bitpack(): data{} {}
        bitpack(std::string bits): data{} {
            size_t rem = size(bits) % bits_per_block;
            if(rem) {
                bits += std::string(bits_per_block - rem, '0');
            }
            for(size_t i = 0, pos = 0; pos < size(bits); i++, pos += bits_per_block) {
                for(size_t j = bits_per_block; j; j--) {
                    data[i] *= 2;
                    data[i] ^= bits[pos + j - 1] == '1';
                }
            }
        }

        bitpack& xor_hint(bitpack const& t, size_t hint) {
            for(size_t i = hint / bits_per_block; i < blocks; i++) {
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

        int operator[](size_t i) const {
            return (data[i / bits_per_block] >> (i % bits_per_block)) & 1;
        }

        std::string to_string() const {
            std::string res(blocks * bits_per_block, '0');
            for(size_t i = 0, pos = 0; i < blocks; i++, pos += bits_per_block) {
                Int block = data[i];
                for(size_t j = 0; j < bits_per_block; j++) {
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
                res += bits_per_block;
                i++;
            }
            if(i < blocks) {
                res += std::countr_zero(data[i]);
            }
            return std::min(res, n);
        }
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_BITPACK_HPP
