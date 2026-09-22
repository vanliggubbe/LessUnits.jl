# LessUnits

[![Build Status](https://github.com/vanliggubbe/LessUnits.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vanliggubbe/LessUnits.jl/actions/workflows/CI.yml?query=branch%3Amain)

A package for converting dimensional quantities into dimensionless ones, supporting
[`Unitful.jl`](https://github.com/JuliaPhysics/Unitful.jl) and
[`DynamicQuantities.jl`](https://github.com/JuliaPhysics/DynamicQuantities.jl).

To use the package, define a tuple of reference quantities that act as the “units” for the corresponding base dimensions you want to work with (e.g. charge, action, frequency, temperature, etc.). You can then convert dimensional quantities to dimensionless values with respect to these references, and
construct the appropriate unit for other dimensional quantities in the same reduced unit system.

## Usage

```julia
import Unitful as UF
using LessUnits

# Reference quantities:
# charge: double elementary charge
# action: reduced Planck constant
# frequency: 1 GHz (with a 2π factor for angular frequency)
# temperature: energy scale (k_B * T)
u = (2UF.u"q", 1UF.u"ħ", 2π * 1UF.u"GHz", 1UF.u"k")

C_ul = unitless(u, 200UF.u"fF")   # unitless capacitance
L_ul = unitless(u, 0.5UF.u"nH")   # unitless inductance
V_ul = unitless(u, 0.2UF.u"mV")   # unitless voltage

Z_uf = sqrt(L_ul / C_ul) * unitof(UF.ElectricalResistance, u)    # unitful impedance
F_uf = inv(2π * sqrt(C_ul * L_ul)) * unitof(UF.u"Hz", u)        # unitful frequency, in Hz

```

### DynamicQuantities

```julia
using LessUnits
import DynamicQuantities as DQ

basis = (2DQ.u"m", 3DQ.u"s")
unitless(basis, 6DQ.u"m")                    # 3.0
unitof(DQ.u"m/s", basis)                    # (2/3) m s⁻¹
unitof(DQ.dimension(DQ.u"m"), basis)         # 2.0 m

lu = LessUnit(basis...)
lu(6DQ.u"m")                               # 3.0
lu(DQ.dimension(DQ.u"m"))                   # 2.0 m

unitless.(basis, [2DQ.u"m", 4DQ.u"m"])      # [1.0, 2.0]
unitless.(basis, DQ.QuantityArray([2DQ.u"m", 4DQ.u"m"]))
```

Broadcasting holds the entire basis fixed, just like
`map(Base.Fix1(unitless, basis), quantities)`, and supports fused expressions and
in-place assignment.

DQ stores dimensions in values, so pass a quantity or dimension value to
`unitof`, rather than a quantity type. Its unit expressions are themselves
quantities: `lu(DQ.u"m")` performs dimensionless conversion; use `unitof` to
construct a reference quantity. The target's magnitude is ignored by `unitof`.
Symbolic DQ inputs (`DQ.us"cm"`, for example) are expanded to base units before
calculation, and dimensional results use expanded dimensions.

### Mixing backends

```julia
import Unitful as UF
import DynamicQuantities as DQ

unitless((2UF.u"m", 3DQ.u"s"), 4UF.u"m/s")  # ≈ 6.0, with a warning
unitof(UF.u"cm", (2DQ.u"m",))                # 2.0 m, as a DQ quantity
```

Mixed inputs select the DQ implementation and convert Unitful quantities using
DQ's Unitful integration. A rate-limited warning reports this conversion.
The integration supports SI physical dimensions and requires Unitful's default
preferred SI units. Unitful-only calls retain Unitful output and requested units.

For dimensional targets, basis dimensions must be linearly independent and span
the target; otherwise an `ArgumentError` is thrown. Dimensionless targets bypass
that check: `unitof` returns numeric unity and `unitless` returns the numerical
dimensionless value.

## Related packages

Similar functional can be found in [`Dimensionless.jl`](https://github.com/martinkosch/Dimensionless.jl).
