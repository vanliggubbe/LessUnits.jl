# Repository guidance

## Development and verification
- Run commands from the repository root:
  - Setup: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
  - Full package tests: `julia --project=. -e 'using Pkg; Pkg.test()'`
  - Direct test-file run after setup: `julia --project=. test/runtests.jl`
- `test/runtests.jl` loads `test/unitful.jl`, `test/dynamicquantities.jl`, and
  `test/mixed.jl`, then checks method ambiguities. There is no configured lint task.
- Focused backend test: `julia --project=. -e 'using LessUnits, Test; import Unitful as UF; import DynamicQuantities as DQ; include("test/dynamicquantities.jl")'`
- Preserve Julia 1.6 compatibility. CI tests Julia 1.6, 1.12, and prerelease
  on Linux x64. DQ 1.1 supports Julia 1.6; DQ 1.2+ requires Julia 1.10.
- `/Manifest*.toml` is ignored. Use `Project.toml` as the authority for
  package version and compatibility, not a local manifest.

## Implementation gotchas
- Actual API order is `unitless(basis, q)` but `unitof(q, basis)`.
- `unitless.(basis, q)` treats the entire basis tuple as a scalar via a
  function-specific `Base.Broadcast.broadcasted` method; keep it lazy for fusion.
- `unitof` and `unitless` in `src/LessUnits.jl` are generated dispatchers: inspect
  the target argument type and `fieldtypes(basis)` with `_isdynamic`, then emit a
  direct backend call. Type-valued targets require unwrapping `Type{T}`.
- Define generation-time classification helpers before the generated methods.
  Keep conversions and logging in runtime backend calls, not the generators.
- `src/unitful.jl` computes dimensional algebra from types using generated
  functions; `src/dynamicquantities.jl` solves it from values at runtime. DQ backend
  identity is type-inferable even though its physical dimensions are not.
  Both use exact rational dimensional algebra. Import the packages as `UF` and `DQ`.
- UF→DQ warnings live in `_to_dynamic(::UF.Quantity)` and
  `_dynamic_dimension(::UF.Dimensions)`, each with `maxlog=1`; suppression is per
  logging site, so one public call can emit two warnings with a fresh logger.
- Normalize Unitful values before DQ conversion to preserve scales of integer
  quantities in non-SI units.
- Use DQ's `dimension`, `ustrip`, `uexpand`, and dimension `keys`/indexing interfaces.
  Call `uexpand` only on symbolic quantities: DQ 1.1 has no ordinary-quantity fallback.
  Compute reference scales separately from dimensions: fractional powers of whole
  DQ quantities can round intermediate exponents with fixed-denominator storage.
- DQ quantity types do not encode physical dimensions; `unitof` needs a value.
  Symbolic DQ inputs are expanded before calculation, including their scale.
- Basis quantities must have independent dimensions, and target dimensions
  must lie in their span; failures raise `ArgumentError`. Matching the
  set of dimension symbols alone is insufficient. Dimensionless targets
  bypass the independence check.
- For Unitful-only inputs, `unitof(UF.u"s", basis)` selects the output units;
  dimension/type overloads construct the reference quantity. Tests check
  the returned unit as well as numerical equivalence.
- `LessUnit` is callable: quantities use `unitless`, while dimension values,
  Unitful units, and type arguments dispatch to `unitof`. DQ unit expressions
  are quantities, so use `unitof` explicitly to construct their reference quantity.
