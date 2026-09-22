module LessUnits

import Unitful as UF
import DynamicQuantities as DQ
import LinearAlgebra

export unitless, unitof, LessUnit

const Basis = Tuple{Vararg{Union{UF.Quantity, DQ.UnionAbstractQuantity}}}

_isunitful(q) = false
_isunitful(::Type{<:UF.Dimensions}) = true
_isunitful(::Type{<:Union{UF.Quantity{T, D, U}, UF.Level{L, S, UF.Quantity{T, D, U}} where {L, S}} where {T, U}}) where {D} = true
_isunitful(::Type{Type{T}}) where {T} = _isunitful(T)

_isdynamic(q) = false
_isdynamic(::Type{<:Union{DQ.UnionAbstractQuantity, DQ.AbstractDimensions}}) = true
_isdynamic(::Type{Type{T}}) where {T} = _isdynamic(T)

"""
    unitof(q, basis::Tuple)

Return the reference quantity (including its scale) for the dimensions specified
by `q`, treating each quantity in `basis` as unity. The magnitude of `q` is ignored.

Unitful targets may be quantities, quantity types, dimensions, or units. For
Unitful-only inputs, a units target such as `UF.u"s"` also selects the output units.
DynamicQuantities targets must be quantity or dimension values: their types do
not encode physical dimensions. Symbolic DQ quantities are expanded to base units.

Mixing Unitful and DynamicQuantities inputs selects the DynamicQuantities backend
and emits a rate-limited warning. Dimensional results then use DQ quantities.

For a dimensional target, throw `ArgumentError` if the basis dimensions are
linearly dependent or the target dimensions cannot be expressed in the basis.
A dimensionless target returns numeric unity without checking basis independence.
"""
@generated function unitof(q, basis :: Basis)
    if _isdynamic(q) || any(_isdynamic, fieldtypes(basis))
        return :(_dq_unitof(q, basis))
    end
    return :(_uf_unitof(q, basis))
end

"""
    unitless(basis::Tuple, q)

Return the dimensionless value of `q`, treating each quantity in `basis` as unity.
The basis and target may use Unitful or DynamicQuantities. Mixed inputs are
converted to DynamicQuantities with a rate-limited warning.

For a dimensional target, throw `ArgumentError` if the basis dimensions are
linearly dependent or the target dimensions cannot be expressed in the basis.
Dimensionless inputs bypass the basis independence check.

Broadcasting treats the entire basis tuple as a scalar: `unitless.(basis, q)`
is equivalent to `map(Base.Fix1(unitless, basis), q)` for arrays and tuples.
For example, with `import DynamicQuantities as DQ`,
`unitless.((2DQ.u"m",), [2DQ.u"m", 6DQ.u"m"])` returns `[1.0, 3.0]`.
"""
@generated function unitless(basis :: Basis, q)
    if _isdynamic(q) || any(_isdynamic, fieldtypes(basis))
        return :(_dq_unitless(basis, q))
    end
    return :(_uf_unitless(basis, q))
end

Base.Broadcast.broadcasted(::typeof(unitless), basis :: Basis, q) =
    Base.Broadcast.broadcasted(unitless, Ref(basis), q)

"""
    LessUnit(basis::Tuple)
    LessUnit(basis...)

Store a reference basis of Unitful and/or DynamicQuantities quantities. Calling
the result on a quantity computes `unitless(basis, q)`; calling it on dimensions,
Unitful units, or a Unitful quantity type computes `unitof(q, basis)`.
DQ unit expressions are quantities, so use `unitof` explicitly to obtain their
reference quantity.
"""
struct LessUnit{T <: Basis}
    basis :: T
end

LessUnit(a :: Any...) = LessUnit(a)

(a :: LessUnit)(b) = unitless(a.basis, b)
(a :: LessUnit)(b :: Union{UF.Dimension, UF.Dimensions, UF.Units, DQ.AbstractDimensions}) = unitof(b, a.basis)
(a :: LessUnit)(b :: Type) = unitof(b, a.basis)

include("unitful.jl")
include("dynamicquantities.jl")

end
