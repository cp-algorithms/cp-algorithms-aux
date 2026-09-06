# Polynomials and series

`poly_t<T>` owns a finite, normalized coefficient vector. It keeps arithmetic,
coefficient access, slicing, shifts by powers of x, reversal and single-point
evaluation. Algorithms are free functions in `cp_algo::math`:

```cpp
#include "cp-algo/math/poly/series.hpp"
using namespace cp_algo::math;
using mint = modint<998244353>;
using polyn = poly_t<mint>;

polyn p({0, 1, 2});
auto q = exp(p, 100);             // p is preserved
p = exp(std::move(p), 100);       // reuse p's storage where possible
```

| Header | Operations |
| --- | --- |
| `base.hpp` | Representation, arithmetic, coefficient operations |
| `inv.hpp` | Inverse modulo x^n |
| `div.hpp` | Division, remainder, `divmod` |
| `calculus.hpp` | `deriv`, `integr` |
| `series.hpp` | Truncated `log`, `exp`, `pow` |
| `sqrt.hpp` | Truncated square root |
| `euclid.hpp` | `gcd`, `inv_mod`, `resultant` |
| `recurrence.hpp` | `min_rec`, `kth_rec`, intervals of inverse coefficients |
| `eval.hpp` | Multipoint `eval`, `inter`, `to_newton` |
| `chirpz.hpp` | Geometric evaluation and interpolation |
| `transform.hpp` | Taylor shift, Borel transforms, correlations, prefix sum |
| `powmod.hpp` | Powers modulo a polynomial or x^m - 1 |
| `compose.hpp` | Composition |

`../poly.hpp` includes all finite-polynomial algorithms. Include `div.hpp` when
using polynomial `/` or `%`. Helpers for Euclidean algorithms and series inversion
are internal; callers generally need `gcd`, `inv_mod`, `min_rec` or `inv`.

The member algorithm API has been removed. Change `p.exp(n)` to `exp(p, n)`,
`p.inv_inplace(n)` to `p = inv(std::move(p), n)`, and `polyn::inter(x, y)` to
`inter(x, y)`. Value arguments accept both copies and moves; read-only inputs use
const references. Coefficient mutators still return references, so return a named
local directly rather than returning `p.mod_xk_inplace(n)` by value.

`n` is the number of coefficients, not the degree. An empty vector represents
zero; truncation to zero coefficients returns zero. Polynomial coefficients must
form a suitable field for division, and calculus requires invertible integer
denominators. FFT-based multiplication uses the existing modular SIMD backend.

## Lazy FPS and Laurent series

`../fps.hpp` supplies `fps<T>`, a shared handle to memoized coefficients. Requests
compute the missing prefix; smaller subsequent requests reuse it. There is no
fixed precision: available memory and the coefficient field bound usable sizes.

```cpp
fps<mint> fibonacci([](size_t n, auto const& prev) {
    return n < 2 ? mint(1) : prev[n - 1] + prev[n - 2];
});
auto product = fibonacci * fibonacci;
auto first = product.prefix(100); // finite poly_t<mint>
auto later = product[10000];      // extend the same caches

auto e = exp(fps<mint>(polyn({0, 1}))); // coefficients of e^x
```

Addition, subtraction, multiplication, division, derivative, integral, logarithm
and exponential are lazy. Modular multiplication, inverse, logarithm and exponential use
relaxed block convolution, taking O(n log² n) arithmetic work through coefficient
n. They retain O(n) coefficients per expression node. Small blocks and non-modint
products use quadratic multiplication. Finite-polynomial algorithms remain the
preferred path when the final precision is known.

Generators receive their previously computed coefficients and must be causal.
Handles share caches and are not thread-safe. Invalid inverse/log/exp constants
and requests that re-enter an unfinished coefficient throw exceptions. Derivative
can request one extra input coefficient. User generators must not outlive objects
they capture by reference.

`../laurent.hpp` represents `x^offset * series`:

```cpp
laurent<mint> p{fps<mint>(polyn({2, 3})), -2}; // 2/x² + 3/x
mint c = p[-1];
auto q = inv(p);                            // offset +2
```

The offset is an explicit lower bound. Inversion requires the coefficient there
to be nonzero; it does not search indefinitely for a leading term. Addition aligns
bounds, multiplication adds them, and `shift(p, k)` multiplies by x^k. `prefix(n)`
returns coefficients starting at the stored offset. Integration rejects a nonzero
x^-1 coefficient because its antiderivative needs a logarithmic term.

## Verification

```sh
oj-verify run -j 1 verify/poly/*.test.cpp
oj-verify run -j 1 verify/fps/*.test.cpp
g++ -std=c++23 -O2 -I. tests/poly.cpp -o /tmp/poly-properties
/tmp/poly-properties
```

This version of `oj-verify` takes files; expand directory globs in the shell.
It skips files whose dependencies have already been verified. A fresh checkout
with an empty timestamp cache forces execution; cached official inputs and
checkers can be reused. Keep before/after builds in separate directories because
verification programs for the same problem share an `a.out` path.

For performance, save each freshly compiled binary as
`<test-directory>-<test-stem>.out` (for example `poly-exp.test.out`). Then alternate
runs on the same cached inputs, with compilation and other verification stopped:

```sh
python tests/bench_poly.py --before /tmp/before --after /tmp/after \
    --repeat 7 --cases 3 --output /tmp/timings.json \
    verify/poly/{exp,inv,log,pow,div}.test.cpp
```

The script pins itself to one available CPU, warms both executables, alternates
order, and records wall time, CPU time, and executable hashes. It uses a blocking
wait so timeout polling does not round short timings up. Timings include
input/output. Keep correctness validation separate: square roots and randomized
algorithms may legitimately produce different correct outputs.
