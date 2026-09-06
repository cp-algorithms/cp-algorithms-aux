#include "cp-algo/math/fft.hpp"
#include <random>
#include <iostream>
using namespace cp_algo::math::fft;

int main() {
    std::mt19937 rng(42);
    auto random_point = [&]() {
        return point(int(rng() % 1001) - 500, int(rng() % 1001) - 500);
    };
    for(size_t n: {4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096, 65536, 1048576}) {
        std::vector<size_t> cuts{0, 1, 2, 3, 4, 5, n / 8, n / 4, n / 2, n / 2 + 1, n - 1, n};
        if(n <= 128) {
            for(size_t k = 0; k <= n; k++) {cuts.push_back(k);}
        }
        std::ranges::sort(cuts);
        cuts.erase(std::unique(begin(cuts), end(cuts)), end(cuts));
        for(size_t z: cuts) {
            z = std::min(z, n);
            auto test = [&]<bool partial>() {
                cvector a(n);
                for(size_t i = 0; i < z; i++) {a.set(i, random_point());}
                auto b = a;
                a.fft<partial>();
                b.fft<partial>(z);
                for(size_t i = 0; i < n; i++) {
                    assert(abs(a.get(i) - b.get(i)) < 1e-7);
                }
                auto inverse = [&]<bool normalize>() {
                    // Arbitrary spectra: inverse pruning cannot assume a zero output tail.
                    cvector c(n);
                    for(size_t i = 0; i < n; i++) {c.set(i, random_point());}
                    auto d = c;
                    c.ifft<partial, normalize>();
                    d.ifft<partial, normalize>(z);
                    double scale = normalize ? 1 : double(partial ? n / flen : n);
                    for(size_t i = 0; i < z; i++) {
                        assert(abs((c.get(i) - d.get(i)) / scale) < 1e-7);
                    }
                };
                inverse.template operator()<true>();
                inverse.template operator()<false>();
            };
            test.template operator()<true>();
            test.template operator()<false>();
        }
    }
    std::cout << "Zero-tail forward and arbitrary-spectrum inverse prefix properties passed\n";
}
