namespace data_structures {
    namespace segment_tree {
        template<typename meta>
        struct segment_tree {
            const int N;
            vector<meta> _meta;

            segment_tree(int n): N(n), _meta(4 * N) {}

            segment_tree(vector<meta> leafs): N(size(leafs)), _meta(4 * N) {
                build(leafs);
            }

            void pull(int v, int l, int r) {
                if(r - l > 1) {
                    _meta[v].pull(_meta[2 * v], _meta[2 * v + 1], l, r);
                }
            }

            void push(int v, int l, int r) {
                if(r - l > 1) {
                    _meta[v].push(&_meta[2 * v], &_meta[2 * v + 1], l, r);
                } else {
                    _meta[v].push(nullptr, nullptr, l, r);
                }
            }

            void build(auto &a, int v, size_t l, size_t r) {
                if(r - l == 1) {
                    if(l < size(a)) {
                        _meta[v] = a[l];
                    }
                } else {
                    size_t m = (l + r) / 2;
                    build(a, 2 * v, l, m);
                    build(a, 2 * v + 1, m, r);
                    pull(v, l, r);
                }
            }

            void build(auto &a) {
                build(a, 1, 0, N);
            }

            void exec_on_segment(int a, int b, auto func, int v, int l, int r) {
                push(v, l, r);
                if(a <= l && r <= b) {
                    func(_meta[v]);
                    push(v, l, r);
                } else if(r <= a || b <= l) {
                    return;
                } else {
                    int m = (l + r) / 2;
                    exec_on_segment(a, b, func, 2 * v, l, m);
                    exec_on_segment(a, b, func, 2 * v + 1, m, r);
                    pull(v, l, r);
                }
            }

            void exec_on_segment(int a, int b, auto func) {
                exec_on_segment(a, b, func, 1, 0, N);
            }
        };

        struct null_meta {
            using meta = null_meta;
            void pull(meta const&, meta const&, int, int) {}
            void push(meta*, meta*, int, int) {}
        };
    }
}
