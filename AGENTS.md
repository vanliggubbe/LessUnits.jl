# Repository guidance

## Development and verification
- Run commands from the repository root:
  - Setup: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
  - Full package tests: `julia --project=. -e 'using Pkg; Pkg.test()'`
  - Direct test-file run after setup: `julia --project=. test/runtests.jl`
- All tests are in one testset in `test/runtests.jl`; there is no custom
  single-test selector or configured lint/format/typecheck task.
- Preserve Julia 1.6 compatibility. CI tests Julia 1.6, 1.12, and prerelease
  on Linux x64; Unitful compatibility is `1`.
- `/Manifest*.toml` is ignored. Use `Project.toml` as the authority for
  package version and compatibility, not a local manifest.

## Implementation gotchas
- Actual API order is `unitless(basis, q)` but `unitof(q, basis)`.
  The `unitless` docstring incorrectly shows the arguments reversed;
  follow the implementation and tests.
- In `src/LessUnits.jl`, `psinv` and the dimensional `unitof` method are
  generated functions. Dimensional algebra uses exact `Rational{Int}`
  arithmetic from type information; emitted expressions apply the
  runtime basis quantities.
- Basis quantities must have independent dimensions, and target dimensions
  must lie in their span; failures raise `ArgumentError`. Matching the
  set of dimension symbols alone is insufficient.
- `unitof(u"s", basis)` converts the result to the requested Unitful units;
  dimension/type overloads construct the reference quantity. Tests check
  the returned unit as well as numerical equivalence.
- `LessUnit` is callable: quantities use `unitless`, while `Dimension`,
  `Units`, and type arguments dispatch to `unitof`.
