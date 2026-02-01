# LessUnits

[![Build Status](https://github.com/vanliggubbe/LessUnits.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vanliggubbe/LessUnits.jl/actions/workflows/CI.yml?query=branch%3Amain)

A package for converting dimensional quantities into dimensionless ones. It is compatible with `Unitful.jl`.

To use the package, define a tuple of reference quantities that act as the “units” for the corresponding base dimensions you want to work with (e.g. charge, action, frequency, temperature, etc.). You can then convert dimensional quantities to dimensionless values with respect to these references, and
construct the appropriate unit for other dimensional quantities in the same reduced unit system.

## Usage

```julia
using Unitful
using LessUnits

# Reference quantities:
# charge: double elementary charge
# action: reduced Planck constant
# frequency: 1 GHz (with a 2π factor for angular frequency)
# temperature: energy scale (k_B * T)
u = (2u"q", 1u"ħ", 2π * 1u"GHz", 1u"k")

C_ul = unitless(200u"fF", u)   # unitless capacitance
L_ul = unitless(0.5u"nH", u)   # unitless inductance
V_ul = unitless(0.2u"mV", u)   # unitless voltage

Z_uf = sqrt(L_ul / C_ul) * unitof(Unitful.ElectricalResistance, u)    # unitful impedance
F_uf = inv(2π * sqrt(C_ul * L_ul)) * unitof(u"Hz", u)                 # unitful frequency, in Hz

```