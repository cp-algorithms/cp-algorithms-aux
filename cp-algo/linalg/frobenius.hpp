#ifndef CP_ALGO_LINALG_FROBENIUS_HPP
#define CP_ALGO_LINALG_FROBENIUS_HPP
#include "matrix.hpp"
#include "../algebra/poly.hpp"
#include <vector>
namespace cp_algo::linalg {    
    template<bool reduce = false>
    auto frobenius_basis(auto const& A) {
        using matrix = std::decay_t<decltype(A)>;
        using base = matrix::base;
        using polyn = algebra::poly_t<base>;
        assert(A.n() == A.m());
        size_t n = A.n();
        struct krylov {
            std::vector<vec<base>> basis; 
            polyn rec;
        };
        std::vector<krylov> blocks;
        std::vector<vec<base>> reduced;
        while(size(reduced) < n) {
            auto generate_block = [&](auto x) {
                krylov block;
                while(true) {
                    vec<base> y = x | vec<base>::ei(n + 1, size(reduced));
                    for(auto &it: reduced) {
                        y.reduce_by(it);
                    }
                    y.normalize();
                    if(vec<base>(y[std::slice(0, n, 1)]) == vec<base>(n)) {
                        block.rec = std::vector<base>(
                            begin(y) + n + size(reduced) - size(block.basis),
                            begin(y) + n + size(reduced) + 1
                        );
                        return std::pair{block, vec<base>(y[std::slice(n, n, 1)])};
                    } else {
                        block.basis.push_back(x);
                        reduced.push_back(y);
                        x = A.apply(x);
                    }
                }
            };
            auto [block, full_rec] = generate_block(vec<base>::random(n));
            if constexpr (reduce) {
                if(vec<base>(full_rec[std::slice(0, size(reduced), 1)]) != vec<base>(size(reduced))) {
                    auto x = block.basis[0];
                    size_t start = 0;
                    for(auto &[basis, rec]: blocks) {
                        polyn cur_rec = std::vector<base>(
                            begin(full_rec) + start, begin(full_rec) + start + rec.deg()
                        );
                        auto shift = cur_rec / block.rec;
                        for(int j = 0; j <= shift.deg(); j++) {
                            x.add_scaled(basis[j], shift[j]);
                        }
                        start += rec.deg();
                    }
                    reduced.erase(begin(reduced) + start, end(reduced));
                    tie(block, full_rec) = generate_block(x.normalize());
                }
            }
            blocks.push_back(block);
        }
        return blocks;
    }
};
#endif // CP_ALGO_LINALG_FROBENIUS_HPP