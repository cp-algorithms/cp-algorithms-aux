Verified on https://judge.yosupo.jp:
- N = 500'000:
  - Convolution, 440ms (https://judge.yosupo.jp/submission/85695)
  - Convolution (mod 1e9+7), 430ms (https://judge.yosupo.jp/submission/85696)
  - Inv of power series, 713ms (https://judge.yosupo.jp/submission/85694)
  - Exp of power series, 2157ms (https://judge.yosupo.jp/submission/85698)
  - Log of power series, 1181ms (https://judge.yosupo.jp/submission/85699)
  - Pow of power series, 3275ms (https://judge.yosupo.jp/submission/85703)
  - Sqrt of power series, 1919ms (https://judge.yosupo.jp/submission/85705)
  - P(x) -> P(x+a), 523ms (https://judge.yosupo.jp/submission/85706)
  - Division of polynomials, 996ms (https://judge.yosupo.jp/submission/85707)
- N = 100'000:
  - Multipoint evaluation, 2161ms (https://judge.yosupo.jp/submission/85709)
  - Polynomial interpolation, 2551ms (https://judge.yosupo.jp/submission/85711)
  - Kth term of Linear Recurrence, 2913ms (https://judge.yosupo.jp/submission/85727)
- N = 50'000:
  - Inv of Polynomials, 1691ms (https://judge.yosupo.jp/submission/85713)
- N = 10'000:
  - Find Linear Recurrence, 346ms (https://judge.yosupo.jp/submission/85025)

---

The main goal of this library is to implement common polynomial functionality in a
reasonable from competitive programming POV complexity, while also doing it in as
straight-forward way as possible.

Therefore, primary purpose of the library is educational and most of constant-time
optimizations that may significantly harm the code readability were not used.

The library is reasonably fast and generally can be used in most problems where
polynomial operations constitute intended solution. However, it is recommended to
seek out other implementations when the time limit is tight or you really want to
squeeze a solution when it is probably not the intended one.
