# Normalize through public interfaces; symbolic quantities must retain their scale.
_to_dynamic(q :: DQ.UnionAbstractQuantity) = q
_to_dynamic(q :: DQ.UnionAbstractQuantity{<:Any, <:DQ.AbstractSymbolicDimensions}) = DQ.uexpand(q)
_to_dynamic(q :: Number) = q
function _to_dynamic(q :: UF.Quantity)
    @warn "Mixing Unitful and DynamicQuantities inputs; converting to DynamicQuantities." maxlog=1
    preferred = UF.upreferred(q)
    # Normalize first: converting an integer-valued centimetre directly can try
    # to store the metre value in an integer. Rational dimensions avoid rounding
    # Unitful exponents to DQ's default fixed denominator.
    T = typeof(UF.ustrip(preferred))
    return convert(DQ.Quantity{T, DQ.Dimensions{Rational{Int}}}, preferred)
end
_to_dynamic(q :: UF.Level) = _to_dynamic(UF.linear(q))

_dynamic_dimension(q :: DQ.UnionAbstractQuantity) = DQ.dimension(_to_dynamic(q))
_dynamic_dimension(d :: DQ.AbstractDimensions) = d
_dynamic_dimension(d :: DQ.AbstractSymbolicDimensions) = DQ.dimension(DQ.uexpand(DQ.Quantity(1.0, d)))
function _dynamic_dimension(d :: UF.Dimensions) 
    @warn "Mixing Unitful and DynamicQuantities inputs; converting to DynamicQuantities." maxlog=1
    convert(DQ.Dimensions{Rational{Int}}, d)
end
_dynamic_dimension(d :: UF.Dimension) = _dynamic_dimension(UF.Dimensions{(d,)}())
_dynamic_dimension(:: Type{D}) where {D <: UF.Dimensions} = _dynamic_dimension(D())
_dynamic_dimension(q :: Union{UF.Quantity, UF.Level, UF.Units}) = _dynamic_dimension(UF.dimension(q))
_dynamic_dimension(:: Type{<:Union{UF.Quantity{T, D, U}, UF.Level{L, S, UF.Quantity{T, D, U}} where {L, S}} where {T, U}}) where {D} = _dynamic_dimension(D)
function _dynamic_dimension(:: Type{<:Union{DQ.UnionAbstractQuantity, DQ.AbstractDimensions}})
    throw(ArgumentError("DynamicQuantities types do not encode physical dimensions; pass a quantity or dimension value"))
end

function _dq_reference(target :: DQ.AbstractDimensions, basis)
    values = map(DQ.ustrip, basis)
    scale = isempty(values) ? 1 : one(promote_type(map(typeof, values)...))
    iszero(target) && return scale
    isempty(basis) && throw(ArgumentError("A dimensional target cannot be expressed in an empty basis"))

    dimensions = map(DQ.dimension, basis)
    names = unique(vcat(collect(keys(target)), map(d -> collect(keys(d)), dimensions)...))
    exponent(d, name) = name in keys(d) ? Rational{Int}(d[name]) : 0//1
    matrix = [exponent(d, name) for name in names, d in dimensions]
    vector = [exponent(target, name) for name in names]

    # Exact normal equations: singularity and span membership are dimensional
    # facts, not floating-point tolerance decisions.
    coefficients = try
        (matrix' * matrix) \ (matrix' * vector)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        throw(ArgumentError("Basis of dimensions `$(dimensions)` is linearly dependent"))
    end
    matrix * coefficients == vector ||
        throw(ArgumentError("Dimensions `$(target)` cannot be expressed in basis `$(dimensions)`"))

    # Raise only the numerical scales. Raising whole DQ quantities can round
    # intermediate exponents with fixed-denominator dimension representations.
    for (value, power) in zip(values, coefficients)
        iszero(power) && continue
        scale *= value ^ power
    end
    return DQ.Quantity(scale, target)
end

_dq_unitof(q, basis) = _dq_reference(_dynamic_dimension(q), map(_to_dynamic, basis))

function _dq_unitless(basis, q)
    quantity = _to_dynamic(q)
    target = DQ.dimension(quantity)
    reference = _dq_reference(target, map(_to_dynamic, basis))
    result = quantity / reference
    iszero(DQ.dimension(result)) || throw(ArgumentError("Conversion did not produce a dimensionless value"))
    return DQ.ustrip(result)
end
