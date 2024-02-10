#ifndef DATA_STRUCTURES_TREAP_HPP
#define DATA_STRUCTURES_TREAP_HPP
#include <chrono>
#include <random>
/* Submissions on Library Judge:
  Range Reverse Range Sum, 558ms - https://judge.yosupo.jp/submission/147860
  Cartesian Tree, 229ms - https://judge.yosupo.jp/submission/147858
  Dynamic Sequence Range Affine Range Sum, 2245ms - https://judge.yosupo.jp/submission/148948
*/
namespace data_structures {
    namespace treap {
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        #define safe(t, op) (t ? t->op : typename remove_reference<decltype(t->op)>::type())
        template<typename meta>
        struct treap_node {

            using node = treap_node;
            using treap = node*;
            meta _meta;
            int prior = rng();
            size_t size = 1;
            treap children[2] = {nullptr, nullptr};
            enum subtree {L, R};

            base check_sum() {
                push();
                return _meta.val + safe(children[L], check_sum()) + safe(children[R], check_sum());
            }

            treap pull() {
                _meta.pull(children[L], children[R]);
                size = 1 + safe(children[L], size) + safe(children[R], size);
                return this;
            }

            treap push() {
                _meta.push(children[L], children[R]);
                return this;
            }

            // set i-th child and pull metadata
            treap set(subtree i, treap t) {
                children[i] = t;
                return pull();
            }

            // push changes and detach the i-th child
            treap cut(subtree i) {
                return children[i];
            }

            static treap merge(treap A, treap B) {
                if(!safe(A, push()) || !safe(B, push())) {
                    return A ? A : B;
                } else if(A->prior < B->prior) {
                    return A->set(R, merge(A->cut(R), B));
                } else {
                    return B->set(L, merge(A, B->cut(L)));
                }
            }

            // return {L, R}, where |L|=k or L=A when |A| < k
            static pair<treap, treap> split(treap A, size_t k) {
                if(!safe(A, push())) {
                    return {nullptr, nullptr};
                } else if(safe(A->children[L], size) >= k) {
                    auto [split_L, split_R] = split(A->cut(L), k);
                    return {split_L, A->set(L, split_R)};
                } else {
                    k -= safe(A->children[L], size) + 1;
                    auto [split_L, split_R] = split(A->cut(R), k);
                    return {A->set(R, split_L), split_R};
                }
            }

            static void exec_on_segment(treap &A, size_t l, size_t r, auto func) {
                auto [LM, R] = split(A, r);
                auto [L, M] = split(LM, l);
                func(M);
                A = merge(L, merge(M, R));
            }

            static void insert(treap &A, size_t pos, treap t) {
                auto [L, R] = split(A, pos);
                A = merge(L, merge(t, R));
            }

            static void erase(treap &A, size_t pos) {
                auto [L, MR] = split(A, pos);
                auto [M, R] = split(MR, 1);
                delete M;
                A = merge(L, R);
            }

            static void exec_on_each(treap &A, auto func) {
                if(A) {
                    exec_on_each(A->children[L], func);
                    func(A);
                    exec_on_each(A->children[R], func);
                }
            }

            treap pull_all() {
                safe(children[L], pull_all());
                safe(children[R], pull_all());
                return pull();
            }

            treap push_all() {
                push();
                safe(children[L], push_all());
                safe(children[R], push_all());
                return this;
            }

            static treap build(vector<treap> nodes) {
                vector<treap> st;
                for(auto cur: nodes) {
                    while(st.size() >= 2 && st[st.size() - 2]->prior > cur->prior) {
                        st.pop_back();
                    }
                    if(!st.empty() && st.back()->prior > cur->prior) {
                        cur->set(L, st.back());
                        st.pop_back();
                    }
                    if(!st.empty() && st.back()->prior < cur->prior) {
                        st.back()->set(R, cur);
                    }
                    st.push_back(cur);
                }
                return st.empty() ? nullptr : st[0]->pull_all();
            }
        };

        struct null_meta {
            void pull(auto const, auto const) {}
            void push(auto&, auto&) {}
        };
    }
}
#endif // DATA_STRUCTURES_TREAP_HPP