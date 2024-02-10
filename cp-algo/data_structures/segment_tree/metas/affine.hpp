#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP
#define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP
#include "base.hpp"
namespace cp_algo::data_structures::segment_tree::metas {
    template<typename base>
    struct affine_meta: base_meta<affine_meta<base>> {
        using meta = affine_meta;
        struct lin {
            base a = 1, b = 0;
            lin() {}
            lin(base a, base b): a(a), b(b){}

            // a * (t.a * x + t.b) + b
            lin operator * (lin const& t) const {
                return lin{a * t.a, a * t.b + b};
            }

            base apply(base x) const {
                return a * x + b;
            }
        };

        base sum = 0;
        lin to_push = {};

        affine_meta() {}
        affine_meta(base sum): sum(sum) {}

        void push(meta *L, meta *R, int l, int r) override {
            if(to_push.a != 1 || to_push.b != 0) {
                sum = to_push.a * sum + to_push.b * (r - l);
                if(r - l > 1) {
                    L->to_push = to_push * L->to_push;
                    R->to_push = to_push * R->to_push;
                }
                to_push = {};
            }
        }

        void pull(meta const& L, meta const& R, int, int) override {
            sum = L.sum + R.sum;
        }
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_AFFINE_HPP