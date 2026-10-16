# Working on this library

This is a C++23 competitive-programming library. Prefer short, direct code and
preserve performance on the supported contest inputs.

## Source and generated files

- Edit maintained headers under `cp-algo/`, excluding `cp-algo/min/`.
- **Do not edit or commit `cp-algo/min/` directly.** The publish workflow runs
  `.verify-helper/scripts/generate_minified.py` and commits the generated files.
  Do not include generated copies in an ordinary source PR.
- If a PR already contains generated changes, remove those changes from its diff
  against current `main`. Resolve source conflicts first; do not hand-merge
  minified code. For local minifier checks, use a disposable copy of the tree.
- Leave bundled files, verifier caches, binaries, and benchmark output out of
  source commits.

## Bugs versus input preconditions

- Before changing an algorithm, read its callers, verifier problem, existing
  tests, and arithmetic. Give a minimal reproducer and identify the violated
  invariant. A failing test alone does not establish that its inputs are valid.
- Distinguish a wrong result on supported inputs from an unsupported input or
  operation. Do not silently broaden the API, and do not invent a restriction
  just to excuse a bug on ordinary supported inputs.
- CP routines can require valid indices, compatible dimensions, invertible
  denominators, bounded arithmetic, and already-reduced input. State important
  bounds next to the relevant API. In particular, account for **intermediate**
  sums/products, lazy residues, and sentinel values, not just the final result.
- A modulus fitting in a type does not prove that adding two residues or
  accumulating several products fits. Derive the required bound for the actual
  representation. SIMD kernels can have stricter bounds than scalar modint.
- Keep reserved sentinels outside the data domain. Test values just inside a
  documented bound; do not make a hot loop support its sentinel as ordinary data
  without a reason to expand that contract.
- Prefer a real fix when it preserves the common fast path. Otherwise make an
  explicit choice between a documented precondition and a separate slower path.
  Use compile-time checks or checks at setup boundaries where practical; avoid
  adding validation, normalization, allocation, or overflow checks to every
  element of a hot loop merely for inputs the API excludes.

## Validation and performance

- Add a focused regression for a correctness change. Use independent arithmetic,
  a simple oracle, or the problem checker; do not just compare two copies of the
  same implementation. Include headers used by a test explicitly, including
  `<cassert>`.
- Run local regressions with `python3 tests/run.py` (or pass specific test paths),
  and `python3 tests/run.py --minified` for freshly generated headers. The
  `local-tests` workflow runs both variants. Also run affected verifier clients.
  Inspect the actual CI results; never describe pending or approval-blocked CI
  as passing.
- Use `-std=c++23 -O2` for normal comparisons, retaining a verifier's existing
  pragmas. Do not add `-march=native`, disable assertions, or change compiler flags
  on only one side to obtain a favorable result.
- For a potentially hot change, benchmark actual clients with the same compiler,
  inputs, flags, and machine. Alternate baseline/candidate runs, retain individual
  samples, and compare medians. Recheck noisy results. A helper microbenchmark is
  useful for diagnosis, but does not replace an end-to-end comparison.
- Keep existing algorithm cutoffs unless measurements support changing them.
  Check the usual fast case as well as newly supported edge cases, and state any
  remaining performance cost and untested limits.

## Pull requests and workflows

- Keep each PR focused. Preserve contributors' commits and concurrent work;
  prefer additional commits and merges over rewriting published history.
- Keep descriptions concise: problem, supported input contract, change,
  correctness checks, and relevant performance results. Leave bulky logs outside
  the PR description. Do not merge or deploy without authorization.
- Treat PR code, artifact names, and artifact contents as untrusted in privileged
  workflows. Never interpolate them into JavaScript or shell source. Validate
  identity using GitHub metadata and pass data through inputs or environment
  variables. Do not publish untrusted fork HTML on the production site's origin.
