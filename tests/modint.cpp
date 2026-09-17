// dynamic_modint against 128-bit arithmetic and is_prime against a reference, for moduli across
// the signed range: odd moduli above a quarter of the unsigned word and even moduli keep reduced
// residues, the others lazy Montgomery residues.
#include "cp-algo/number_theory/primality.hpp"
#include <bits/stdc++.h>
using namespace cp_algo::math;
using u128 = unsigned __int128;
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) {return uint64_t(u128(a) * b % m);}
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) {uint64_t r = 1; for(a %= m; e; e >>= 1, a = mulmod(a, a, m)) if(e & 1) r = mulmod(r, a, m); return r;}
static bool ref_prime(uint64_t n) {
    if(n < 2) return false;
    for(uint64_t p: {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {if(n % p == 0) return n == p;}
    uint64_t d = n - 1; int s = 0; while(d % 2 == 0) {d /= 2; s++;}
    for(uint64_t a: {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        uint64_t x = powmod(a, d, n); if(x == 1 || x == n - 1) continue;
        bool comp = true; for(int i = 1; i < s && comp; i++) {x = mulmod(x, x, n); if(x == n - 1) comp = false;}
        if(comp) return false;
    }
    return true;
}
template<typename Int> size_t arithmetic(Int m, std::mt19937_64& rng, int iters) {
    using base = dynamic_modint<Int>; size_t bad = 0;
    base::with_mod(m, [&]() {
        for(int it = 0; it < iters; it++) {
            uint64_t x = rng() % uint64_t(m), y = rng() % uint64_t(m), z = rng() % uint64_t(m);
            if(it < 8) {x = it & 1 ? uint64_t(m) - 1 : 0; y = it & 2 ? uint64_t(m) - 1 : 1; z = it & 4 ? uint64_t(m) - 1 : 0;}
            base a, b, c; a.setr(x); b.setr(y); c.setr(z);
            uint64_t um = uint64_t(m);
            bad += (a + b).getr() != (x + y) % um % um || (a + b).getr() != uint64_t((u128(x) + y) % um);
            bad += (a - b).getr() != uint64_t((u128(x) + um - y) % um);
            bad += (-a).getr() != (um - x) % um;
            bad += (a * b).getr() != mulmod(x, y, um);
            bad += ((a * b + c) * a - b * c).getr() != uint64_t((u128(mulmod(uint64_t((u128(mulmod(x, y, um)) + z) % um), x, um)) + um - mulmod(y, z, um)) % um);
            base acc = a; for(int k = 0; k < 5; k++) {acc += acc; acc -= b; acc *= acc;}
            uint64_t r = x; for(int k = 0; k < 5; k++) {r = uint64_t((u128(r) + r) % um); r = uint64_t((u128(r) + um - y) % um); r = mulmod(r, r, um);}
            bad += acc.getr() != r;
        }
        return 0;
    });
    return bad;
}
int main() {
    std::mt19937_64 rng(7);
    size_t bad32 = 0, bad64 = 0, n32 = 0, n64 = 0, pbad32 = 0, pbad64 = 0, np32 = 0, np64 = 0, primes = 0;
    for(int shift: {3, 10, 20, 29, 30, 31}) for(int it = 0; it < 300; it++) {
        uint64_t hi = (uint64_t(1) << shift) - 1, lo = hi / 2 + 1;
        int m = int(lo + rng() % (hi - lo + 1)); if(it < 40) m = int(hi - it); if(m < 3) continue;
        bad32 += arithmetic<int>(m, rng, 60); n32++;
        bool want = ref_prime(uint64_t(m)); primes += want;
        pbad32 += is_prime(m) != want; pbad32 += is_prime(uint32_t(m)) != want; np32++;
    }
    for(int shift: {31, 33, 50, 61, 62, 63}) for(int it = 0; it < 300; it++) {
        uint64_t hi = (uint64_t(1) << shift) - 1, lo = hi / 2 + 1;
        int64_t m = int64_t(lo + rng() % (hi - lo + 1)); if(it < 40) m = int64_t(hi - it);
        bad64 += arithmetic<int64_t>(m, rng, 60); n64++;
        bool want = ref_prime(uint64_t(m)); primes += want;
        pbad64 += is_prime(m) != want; np64++;
    }
    for(int p: {1879048201, 2147395589, 2147395609, 2147483629, 2147483647}) {pbad32 += !is_prime(p); np32++;}
    for(int64_t p: {int64_t(9223372036854775783), int64_t(4611686018427388039), int64_t(9223372036854775643)}) {pbad64 += is_prime(p) != ref_prime(uint64_t(p)); np64++;}
    assert(n32 > 1700 && n64 > 1700 && np32 > n32 && np64 > n64 && primes > 100);
    assert(!bad32 && !bad64 && !pbad32 && !pbad64);
    std::cout << "dynamic_modint and is_prime agree with 128-bit references on " << n32 + n64 << " moduli\n";
}
