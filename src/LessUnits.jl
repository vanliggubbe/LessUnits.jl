module LessUnits

import Unitful: Dimensions, Dimension, Quantity, @u_str, dimension, uconvert, NoUnits
import LinearAlgebra: I

export unitless

unravel(a :: Dimension{T}) where {T} = (T => a.power)
unravel(:: Dimensions{T}) where {T} = map(unravel, T)
unravel(:: Type{Dimensions{T}}) where {T} = map(unravel, T)

@generated function dim_mat(basis :: Dimensions...)
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
    return quote
        $(index), $(reduce(hcat, vecs))
    end
end

"""
    unitless(q :: Quantity, basis :: Tuple{Vararg{Quantity}})

Returns dimensionless value of quantity `q`, assuming each element of `basis` to be normalized to unity.
Throws `ArgumentError` if elements of `basis` are not independent or `q` cannot be expressed through `basis` units.
"""
@generated function unitless(q :: Quantity, basis :: Tuple{Vararg{Quantity}})
    index, basis_dim = dim_mat(dimension.(fieldtypes(basis))...)
    #index, basis_dim = dim_mat(dimension.(basis)...)
    q_dim = zeros(Rational{Int}, size(basis_dim)[1])
    for (dim, pow) in unravel(dimension(q))
        if !haskey(index, dim)
            return :(throw(ArgumentError("Quantity of dimension `$(dimension(q))` cannot be expanded over `$(dimension.(basis))`")))
        end
        q_dim[index[dim]] = pow
    end
    # solve system to expand quantity over basis
    A = basis_dim' * basis_dim
    b = basis_dim' * q_dim
    n = length(b)
    # Gaussian elimination
    for i in 1 : (n - 1)
        j = i + argmax(map(abs, A[(i + 1) : end, i]))
        A[i, :], A[j, :] = A[j, :], A[i, :]
        b[i], b[j] = b[j], b[i]
        if iszero(A[i, i])
            return :(throw(ArgumentError("Basis of dimensions `$(dimension.(basis))` is linearly dependent")))
        end
        for k in (i + 1) : n
            tmp = A[k, i] / A[i, i]
            A[k, i : end] -= tmp * A[i, i : end]
            b[k] -= tmp * b[i]
        end
    end
    x = zero(b)
    # solving right-triangular
    for i in n : -1 : 1
        x[i] = b[i] / A[i, i]
        b[1 : (i - 1)] -= A[1 : (i - 1), i] * x[i]
    end

    if basis_dim * x != q_dim
        return :(throw(ArgumentError("Quantity of dimension `$(dimension(q))` cannot be expanded over `$(dimension.(basis))`")))
    end
    ret = :(q)
    for (i, pow) in enumerate(x)
        if !iszero(pow)
            ret = :($(ret) * (basis[$(i)] ^ ($(-pow))))
        end
    end
    return :(uconvert(NoUnits, $(ret)))
end

end
