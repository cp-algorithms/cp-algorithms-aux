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

            void exec_on_segment(int a, int b, auto func, auto proceed, auto stop, int v, int l, int r) {
                push(v, l, r);
                if(r <= a || b <= l || stop(_meta[v])) {
                    return;
                } else if(a <= l && r <= b && proceed(_meta[v])) {
                    func(_meta[v]);
                    push(v, l, r);
                } else {
                    int m = (l + r) / 2;
                    exec_on_segment(a, b, func, proceed, stop, 2 * v, l, m);
                    exec_on_segment(a, b, func, proceed, stop, 2 * v + 1, m, r);
                    pull(v, l, r);
                }
            }

            static constexpr auto default_true = [](auto const&){return true;};
            static constexpr auto default_false = [](auto const&){return false;};

            void exec_on_segment(int a, int b, auto func, auto proceed, auto stop) {
                exec_on_segment(a, b, func, proceed, stop, 1, 0, N);
            }

            void exec_on_segment(int a, int b, auto func) {
                exec_on_segment(a, b, func, default_true, default_false);
            }
        };

        namespace metas {
            struct null_meta {
                using meta = null_meta;
                void pull(meta const&, meta const&, int, int) {}
                void push(meta*, meta*, int, int) {}
            };

            struct chmin_chmax_add_meta {
                static constexpr int64_t inf = 1e12;

                using meta = chmin_chmax_sum_meta;
                int64_t sum = 0, add = 0;

                template<typename Comp>
                struct data {
                    int64_t val;
                    int64_t count = 1;
                    int64_t second = max(inf, -inf, comp);
                    static const Comp comp;

                    data combine(data const& t) const {
                        return comp(val, t.val) ? data{val, count, min(second, t.val, comp)}
                             : comp(t.val, val) ? data{t.val, t.count, min(t.second, val, comp)}
                             : data{val, count + t.count, min(second, t.second, comp)};
                    }

                    void add(int64_t b) {
                        val += b;
                        second += b;
                    }

                    int64_t normalize(int64_t L, int64_t R) {
                        int64_t old_val = val;
                        val = clamp(val, L, R);
                        second = clamp(second, L, R);
                        return count * (val - old_val);
                    }

                    bool stop(int64_t b) const {
                        return !comp(val, b);
                    }
                    bool proceed(int64_t b) const {
                        return comp(b, second);
                    }
                };
                data<less<int64_t>> mn = {sum};
                data<greater<int64_t>> mx = {sum};
                int64_t chmin = inf, chmax = -inf;

                void pull(meta const& L, meta const& R, int, int) {
                    sum = L.sum + R.sum;
                    mn = L.mn.combine(R.mn);
                    mx = L.mx.combine(R.mx);
                }

                void push(meta &t) {
                    t.add += add; t.chmin += add; t.chmax += add;
                    t.chmin = clamp(t.chmin, chmax, chmin);
                    t.chmax = clamp(t.chmax, chmax, chmin);
                }

                void push(meta* L, meta* R, int l, int r) {
                    if(r - l > 1) {
                        push(*L);
                        push(*R);
                    }
                    if(add) {
                        sum += (r - l) * add;
                        mn.add(add);
                        mx.add(add);
                    }
                    bool same = mn.val == mx.val;
                    auto to_add = mn.normalize(chmax, chmin) + mx.normalize(chmax, chmin);
                    sum += same ? to_add / 2 : to_add;
                    if(mn.val == mx.val) {
                        mx = {mx.val, r - l};
                        mn = {mn.val, r - l};
                    }
                    add = 0;
                    chmin = inf;
                    chmax = -inf;
                }

                static auto proceed_chmin(int64_t b) {
                    return [b](meta const& t) {return t.mx.proceed(b);};
                }
                static auto stop_chmin(int64_t b) {
                    return [b](meta const& t) {return t.mx.stop(b);};
                }
                static auto proceed_chmax(int64_t b) {
                    return [b](meta const& t) {return t.mn.proceed(b);};
                }
                static auto stop_chmax(int64_t b) {
                    return [b](meta const& t) {return t.mn.stop(b);};
                }
            };
        }
    }
}
