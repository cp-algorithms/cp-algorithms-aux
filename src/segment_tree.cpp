#include <bits/stdc++.h>

#define int int64_t

using namespace std;

template<int m>
struct modular {
    int r;
    constexpr modular(): r(0) {}
    constexpr modular(int64_t rr): r(rr % m) {if(r < 0) r += m;}
    modular inv() const {return bpow(*this, m - 2);}
    modular operator - () const {return r ? m - r : 0;}
    modular operator * (const modular &t) const {return (int64_t)r * t.r % m;}
    modular operator / (const modular &t) const {return *this * t.inv();}
    modular operator += (const modular &t) {r += t.r; if(r >= m) r -= m; return *this;}
    modular operator -= (const modular &t) {r -= t.r; if(r < 0) r += m; return *this;}
    modular operator + (const modular &t) const {return modular(*this) += t;}
    modular operator - (const modular &t) const {return modular(*this) -= t;}
    modular operator *= (const modular &t) {return *this = *this * t;}
    modular operator /= (const modular &t) {return *this = *this / t;}
    
    bool operator == (const modular &t) const {return r == t.r;}
    bool operator != (const modular &t) const {return r != t.r;}
    
    explicit operator int() const {return r;}
    int64_t rem() const {return 2 * r > m ? r - m : r;}
};

template<int T>
istream& operator >> (istream &in, modular<T> &x) {
    return in >> x.r;
}

template<int T>
ostream& operator << (ostream &out, modular<T> const& x) {
    return out << x.r;
}

const int mod = 998244353;
using base = modular<mod>;

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

        struct affine_meta: null_meta {
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

            void push(meta *L, meta *R, int l, int r) {
                if(to_push.a != 1 || to_push.b != 0) {
                    sum = to_push.a * sum + to_push.b * (r - l);
                    if(r - l > 1) {
                        L->to_push = to_push * L->to_push;
                        R->to_push = to_push * R->to_push;
                    }
                    to_push = {};
                }
            }

            void pull(meta const& L, meta const& R, int, int) {
                sum = L.sum + R.sum;
            }
        };
    }
}

using namespace data_structures::segment_tree;

void solve() {
    int n, q;
    cin >> n >> q;
    vector<affine_meta> a(n);
    for(int i = 0; i < n; i++) {
        int ai;
        cin >> ai;
        a[i] = {ai};
    }
    segment_tree<affine_meta> me(a);
    while(q--) {
        int t;
        cin >> t;
        if(t == 0) {
            int l, r, b, c;
            cin >> l >> r >> b >> c;
            me.exec_on_segment(l, r, [&](auto& meta) {
                meta.to_push = affine_meta::lin(b, c) * meta.to_push;
            });
        } else {
            int l, r;
            cin >> l >> r;
            base ans = 0;
            me.exec_on_segment(l, r, [&](auto meta) {
                ans += meta.sum;
            });
            cout << ans << "\n";
        }
    }
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t; 
    t = 1;// cin >> t;
    while(t--) {
        solve();
    }
}
 
