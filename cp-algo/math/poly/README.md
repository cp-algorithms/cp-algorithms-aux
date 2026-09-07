# Polynomials and series

`poly_t<T>` owns a finite, normalized coefficient vector. It keeps arithmetic,
coefficient access, slicing, shifts by powers of x, reversal and single-point
evaluation. Algorithms are free functions in `cp_algo::math`:

```cpp
#include "cp-algo/math/poly/series/exp.hpp"
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
| `series/inv.hpp` | Inverse modulo x^n |
| `div.hpp` | Division, remainder, `divmod` |
| `calculus.hpp` | `deriv`, `integr` |
| `series/{log,exp,pow}.hpp` | Truncated `log`, `exp`, `pow`; `series.hpp` includes all three |
| `series/sqrt.hpp` | Truncated square root |
| `sparse/{inv,log,exp,pow,sqrt}.hpp` | Sparse algorithms; `sparse.hpp` includes all five |
| `euclid.hpp` | `gcd`, `inv_mod`, `resultant` |
| `recurrence.hpp` | `min_rec`, `kth_rec`, intervals of inverse coefficients |
| `eval.hpp` | Multipoint `eval`, `inter`, `to_newton` |
| `chirpz.hpp` | Geometric evaluation and interpolation |
| `transform.hpp` | Taylor shift, Borel transforms, correlations, prefix sum |
| `powmod.hpp` | Powers modulo a polynomial or x^m - 1 |
| `compose.hpp` | Composition |

The `series/` directory groups finite power-series operations, and `sparse/` groups
their sparse counterparts. The root holds representation, polynomial arithmetic,
evaluation, recurrence and transformation algorithms.

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

Sparse algorithms accept the same `poly_t<T>` and return a dense truncated result.
They collect the nonzero input coefficients and use O(nk) recurrences for k nonzero
terms. Use them for sparse coefficients, even when the degree is large. As with
finite-polynomial calculus, set `CP_ALGO_MAXN` above the requested precision for
the cached integer inverses. These tables require a fixed modulus at least as
large as `CP_ALGO_MAXN`, with the requested precision below the characteristic.

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
auto filtered = fibonacci * polyn({1, 2, 3}); // one known input

auto e = exp(fps<mint>(polyn({0, 1}))); // coefficients of e^x
```

Addition, subtraction, multiplication, division, derivative, integral, logarithm
and exponential are lazy. `pow(p, k)` supports nonnegative integer powers, including
leading zeros and zero generators. It discovers only as many leading coefficients
as the requested output needs; exponent zero never evaluates the input. Power
uses `p q' = k p' q`, retaining the product strategy for known polynomial inputs.
 Modular multiplication, inverse, logarithm and exponential use
relaxed block convolution, taking O(n log² n) arithmetic work through coefficient
n. They retain O(n) coefficients per expression node. Small blocks and non-modint
products use quadratic multiplication. FPS constructed from a polynomial remember
that their input is fully known. Products with at most 16 nonzero positive-degree
terms use a sparse recurrence. Other known inputs use semi-relaxed convolution,
caching their FFT prefixes and computing only the needed middle products.
The same path is used by `inv`, `log`, `exp`, and `pow` when their input is known. Two
generator-backed inputs still use fully relaxed multiplication. Neither convolution
path reads unknown coefficients ahead.
Known polynomial caches stay immutable; queries beyond their degree return zero.
Lazy calculus reuses cached integer inverses within `CP_ALGO_MAXN` for fixed
moduli large enough for the table. It keeps ordinary division for dynamic/small
moduli and beyond the table so requests can continue growing.
Finite-polynomial algorithms remain the preferred path when the final precision is known.

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
oj-verify run -j 1 $(rg --files verify/poly -g '*.test.cpp')
oj-verify run -j 1 verify/fps/*.test.cpp
g++ -std=c++23 -O2 -I. tests/poly.cpp -o /tmp/poly-properties
/tmp/poly-properties
g++ -std=c++23 -O2 -I. tests/fps.cpp -o /tmp/fps-properties
/tmp/fps-properties
g++ -std=c++23 -O2 -I. tests/poly_sparse.cpp -o /tmp/poly-sparse-properties
/tmp/poly-sparse-properties
g++ -std=c++23 -O2 -I. tests/poly_series.cpp -o /tmp/poly-series-properties
/tmp/poly-series-properties
g++ -std=c++23 -O2 -I. tests/power.cpp -o /tmp/power-properties
/tmp/power-properties
```

This version of `oj-verify` takes files. The `rg` command includes the sparse
verifiers in `verify/poly/sparse/`. Expand other directory globs in the shell.
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

## Finite-series algorithms

Dense exponential and square root maintain the result and its reciprocal across
Newton doublings, using only the necessary middle and low products. The final
step is truncated to the requested precision and skips the unused reciprocal
update. Logarithm solves `p*q = p'` in two halves with a half-length inverse.
Retained transforms belong to unchanged coefficient buffers; reduced product
coefficients are transformed again before use in another modular product.

Interpolation evaluates the derivative on the existing product tree and uses
`bulk_invs` for its scalar weights. Its remainder-tree algorithm is unchanged.
Polynomial division extracts only the leading reversed slice and recovers only
coefficients below the divisor's degree in the remainder. Recurrence queries also
truncate their final inverse and product to the requested coefficient.

## Powering and squaring

`powmod` and `powmod_circular` preserve the FFT squaring shortcut when the
operation squares one operand. Small self-products also avoid a temporary
copy and combine symmetric coefficient pairs. Short wrapped-tail products and
recursive large convolutions also reuse a square's direct transforms and
modular splits. Polynomial modular powers
remain binary: windowing did not consistently improve their measured consumers.

The dense matrix `pow` method uses three-bit windows, precomputing odd powers
when the exponent's bit pattern makes this cheaper than binary powering.
This mechanism is also available explicitly:

```cpp
auto result = bpow<3>(value, exponent, identity, operation);
```

The overloads with an explicit identity accept a compile-time window from 1 to 6.
Ordinary `bpow` remains binary by default; windowing needs extra stored values
and does not uniformly help inexpensive scalar multiplication. The two-argument
`bpow(value, exponent)` retains its existing interface and binary algorithm.
Both algorithms require a nonnegative exponent fitting in 64 bits.
