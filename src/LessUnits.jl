module LessUnits

import Unitful: Dimensions, Dimension, Quantity, @u_str, dimension, uconvert, NoUnits, Level
import LinearAlgebra: I

export unitless, unitof

unravel(a :: Dimension{T}) where {T} = (T => a.power)
unravel(:: Dimensions{T}) where {T} = map(unravel, T)
unravel(:: Type{Dimensions{T}}) where {T} = map(unravel, T)

@generated function psinv(basis :: Dimensions...)
    index = Dict{Symbol, Int}()
    vecs = Vector{Rational{Int}}[]
    # expand over SI fundamental units
    for el in basis
        vec = Rational{Int}[]
        for (dim, pow) in unravel(el)
            if !haskey(index, dim)
                index[dim] = length(index) + 1
            end
            idx = index[dim]
            if length(vec) < idx
                append!(vec, zeros(Rational{Int}, idx - length(vec)))
            end
            vec[idx] = pow
        end
        push!(vecs, vec)
    end
    for vec in vecs
        append!(vec, zeros(Rational{Int}, length(index) - length(vec)))
    end
    D = reduce(hcat, vecs)
    
    # find pseudoinverse
    A = D' * D
    B = Matrix(D')
    n = size(A)[1]
    # Gaussian elimination
    for i in 1 : (n - 1)
        j = i + argmax(map(abs, A[(i + 1) : end, i]))
        A[i, :], A[j, :] = A[j, :], A[i, :]
        B[i, :], B[j, :] = B[j, :], B[i, :]
        if iszero(A[i, i])
            return :(throw(ArgumentError("Basis of dimensions `$(basis)` is linearly dependent")))
        end
        for k in (i + 1) : n
            tmp = A[k, i] / A[i, i]
            A[k, i : end] -= tmp * A[i, i : end]
            B[k, :] -= tmp * B[i, :]
        end
    end
    # solving right-triangular
    for i in n : -1 : 1
        if iszero(A[i, i])
            return :(throw(ArgumentError("Basis of dimensions `$(basis)` is linearly dependent")))
        end
        B[i, :] /= A[i, i]
        B[1 : (i - 1), :] -= A[1 : (i - 1), i] * transpose(B[i, :])
    end

    return quote
        $(index), $(D), $(B)
    end
end

_number(:: Quantity{T, D, U}) where {T, D, U} = T
_number(:: Type{Quantity{T, D, U}}) where {T, D, U} = T

"""
    unitof(q, basis :: Tuple{Vararg{Quantity}})

Returns unit of dimensions specified by `q`, assuming each element of `basis` corresponds to unity.
Throws `ArgumentError` if elements of `basis` are not independent or `q` cannot be expressed through `basis` units.
"""
function unitof end

unitof(:: Dimensions{()}, basis :: Tuple{Vararg{Quantity}}) = one(promote_type(_number.(typeof.(basis))...))

@generated function unitof(q :: Dimensions, basis :: Tuple{Vararg{Quantity}})
    index, basis_dim, basis_inv = psinv(dimension.(fieldtypes(basis))...)
    q_dim = zeros(Rational{Int}, size(basis_dim)[1])
    for (dim, pow) in unravel(q)
        if !haskey(index, dim)
            return :(throw(ArgumentError("Quantity of dimension `$(q)` cannot be expanded over `$(dimension.(basis))`")))
        end
        q_dim[index[dim]] = pow
    end
    x = basis_inv * q_dim
    if basis_dim * x != q_dim
        return :(throw(ArgumentError("Quantity of dimension `$(q)` cannot be expanded over `$(dimension.(basis))`")))
    end
    ret = :($(one(promote_type(_number.(fieldtypes(basis))...))))
    for (i, pow) in enumerate(x[1 : end])
        if !iszero(pow)
            ret = :($(ret) * (basis[$(i)] ^ ($(pow))))
        end
    end
    return ret
end

unitof(
    :: Type{<: Union{Quantity{T, D, U}, Level{L, S, Quantity{T, D, U}} where {L, S}} where {T, U}},
    basis :: Tuple{Vararg{Quantity}}
) where {D} = unitof(D, basis)

unitof(:: T, basis :: Tuple{Vararg{Quantity}}) where {T <: Union{Quantity, Level}} = unitof(T, basis)

"""
    unitless(q, basis :: Tuple{Vararg{Quantity}})

Returns dimensionless value of quantity `q`, assuming each element of `basis` to be normalized to unity.
Throws `ArgumentError` if elements of `basis` are not independent or `q` cannot be expressed through `basis` units.
"""
unitless(q, basis :: Tuple{Vararg{Quantity}}) = uconvert(NoUnits, q / unitof(dimension(q), basis))    

end
