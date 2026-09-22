unravel(a :: UF.Dimension{T}) where {T} = (T => a.power)
unravel(:: UF.Dimensions{T}) where {T} = map(unravel, T)
unravel(:: Type{UF.Dimensions{T}}) where {T} = map(unravel, T)

@generated function psinv(basis :: UF.Dimensions...)
    index = Dict{Symbol, Int}()
    vecs = Vector{Rational{Int}}[]
    # Expand over fundamental dimensions.
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
    D = isempty(vecs) ? zeros(Rational{Int}, 0, 0) : reduce(hcat, vecs)

    # Find pseudoinverse.
    A = D' * D
    B = Matrix(D')
    n = size(A)[1]
    # Gaussian elimination.
    for i in 1 : (n - 1)
        j = i - 1 + argmax(map(abs, A[i : end, i]))
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
    # Solve right-triangular system.
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

_number(:: UF.Quantity{T, D, U}) where {T, D, U} = T
_number(:: Type{UF.Quantity{T, D, U}}) where {T, D, U} = T

_uf_unitof(:: UF.Dimensions{()}, basis :: Tuple{Vararg{UF.Quantity}}) =
    isempty(basis) ? 1 : one(promote_type(_number.(typeof.(basis))...))

@generated function _uf_unitof(q :: UF.Dimensions, basis :: Tuple{Vararg{UF.Quantity}})
    index, basis_dim, basis_inv = psinv(UF.dimension.(fieldtypes(basis))...)
    q_dim = zeros(Rational{Int}, size(basis_dim)[1])
    for (dim, pow) in unravel(q)
        if !haskey(index, dim)
            return :(throw(ArgumentError("Quantity of dimension `$(q)` cannot be expanded over `$(UF.dimension.(basis))`")))
        end
        q_dim[index[dim]] = pow
    end
    x = basis_inv * q_dim
    if basis_dim * x != q_dim
        return :(throw(ArgumentError("Quantity of dimension `$(q)` cannot be expanded over `$(UF.dimension.(basis))`")))
    end
    ret = :($(one(promote_type(_number.(fieldtypes(basis))...))))
    for (i, pow) in enumerate(x[1 : end])
        if !iszero(pow)
            ret = :($(ret) * (basis[$(i)] ^ ($(pow))))
        end
    end
    return ret
end

_uf_unitof(u :: UF.Units{N, D, A}, basis) where {N, D, A} = UF.uconvert(u, _uf_unitof(D, basis))

_uf_unitof(
    :: Type{<: Union{UF.Quantity{T, D, U}, UF.Level{L, S, UF.Quantity{T, D, U}} where {L, S}} where {T, U}},
    basis
) where {D} = _uf_unitof(D, basis)

_uf_unitof(:: T, basis) where {T <: Union{UF.Quantity, UF.Level}} = _uf_unitof(T, basis)
_uf_unitof(d :: UF.Dimension, basis) = _uf_unitof(UF.Dimensions{(d,)}(), basis)
_uf_unitof(:: Type{D}, basis) where {D <: UF.Dimensions} = _uf_unitof(D(), basis)

_uf_unitless(basis, q) = UF.uconvert(UF.NoUnits, q / _uf_unitof(UF.dimension(q), basis))
