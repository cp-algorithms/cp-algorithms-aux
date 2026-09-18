// The products must not depend on the caller's allocator zero-initializing what a vector grows
// into: every coefficient a result needs is written, including the padding past its true length.
#include "cp-algo/math/fft.hpp"
#include <random>
#include <cstring>
#include <cstdio>
#include <cassert>
using namespace cp_algo;
using namespace cp_algo::math;
// Grows without initializing, and fills fresh blocks with a pattern, so any reliance on
// zero-initialized growth shows up as a wrong answer.
template<class T> struct poison_alloc: big_alloc<T> {
    using big_alloc<T>::big_alloc;
    poison_alloc() = default;
    template<class U> poison_alloc(poison_alloc<U> const&) noexcept {}
    template<class U> struct rebind {using other = poison_alloc<U>;};
    [[nodiscard]] T* allocate(std::size_t n) {
        T* p = big_alloc<T>::allocate(n);
        std::memset(static_cast<void*>(p), 0xCC, n * sizeof(T));
        return p;
    }
    template<class U, class... Args> void construct(U* p, Args&&... args) {
        if constexpr(sizeof...(Args) == 0) {::new(static_cast<void*>(p)) U;}
        else {::new(static_cast<void*>(p)) U(std::forward<Args>(args)...);}
    }
};
template<class T> using poison_vector = std::vector<T, poison_alloc<T>>;
using base = modint<998244353>;
int bad = 0;
template<class V> V make(size_t n, std::mt19937_64& g) {
    V v; v.resize(n);
    for(auto &x: v) {x = base(g() % base::mod());}
    return v;
}
void check(const char* what, size_t as, size_t bs, size_t k, bool truncate, std::mt19937_64& g) {
    auto a = make<big_vector<base>>(as, g), b = make<big_vector<base>>(bs, g);
    size_t need = truncate ? k : as + bs - 1;
    big_vector<base> want(need);
    for(size_t i = 0; i < as; i++) {for(size_t j = 0; j < bs && i + j < need; j++) {want[i + j] += a[i] * b[j];}}
    poison_vector<base> pa, pb;
    pa.resize(as); pb.resize(bs);
    std::copy(a.begin(), a.end(), pa.begin());
    std::copy(b.begin(), b.end(), pb.begin());
    if(truncate) {fft::mul_truncate(pa, pb, k);} else {fft::mul(pa, pb);}
    bool ok = pa.size() == need;
    for(size_t i = 0; ok && i < need; i++) {ok = pa[i] == want[i];}
    if(!ok) {printf("  relies on zeros: %-14s as=%-6zu bs=%-6zu k=%-6zu\n", what, as, bs, k); bad++;}
}
int main() {
    std::mt19937_64 g(5);
    for(auto [as, bs]: {std::pair<size_t,size_t>{3,4}, {10,10}, {63,63}, {64,64}, {100,70},
                        {1000,1000}, {5000,3}, {70000,100}, {1<<16,1<<16}}) {
        check("mul", as, bs, 0, false, g);
        for(size_t k: {size_t(1), (as + bs) / 2, as + bs - 1, as + bs + 50}) {
            check("mul_truncate", as, bs, k, true, g);
        }
    }
    // the overlap-add path for a long operand against a short one
    check("unbalanced", 1 << 20, 100, 0, false, g);
    assert(!bad);
    printf("No product relies on zero-initialized growth\n");
}
