@def title = "Studynote on Single-Matrix Face Detection and Exact Kernel Recovery for PSD Matrices in  ClusteredLowRankSolver.jl."   
@def published = "February 15, 2026"   
@def tags =["programming", "Julia"]  

# Studynote on Single-Matrix Face Detection and Exact Kernel Recovery for PSD Matrices in  ClusteredLowRankSolver.jl

**Shuvomoy Das Gupta**
*February 15, 2026"*

**Table of contents**
\toc

## Abstract

In this study note, we will discuss how to use the rounding procedure described in `ClusteredLowRankSolver.jl`. To that goal, we will use a simplified and self-contained version of the rounding procedure described in the file `https://github.com/nanleij/ClusteredLowRankSolver.jl/blob/main/src/rounding.jl`. 

The simplified code is provided in Section 0. Please note that this is nothing but a simplified version of 

`https://github.com/nanleij/ClusteredLowRankSolver.jl/blob/main/src/rounding.jl`

and the modification that we make here are for the sake of illustration, for complete implementation  see `ClusteredLowRankSolver.jl`. The simplified code implements a single-matrix rounding/face-detection pipeline for a numerical symmetric PSD matrix $\widehat{Y}\in\mathbb{S}^n$. 

We also try to explain the mathematics behind the single-matrix rounding / face-detection pipeline implemented in `rounding_modified.jl`. The pipeline takes as input a numerical symmetric PSD matrix $\widehat{Y}\in\mathbb{S}^n$ and outputs a matrix $Q\in\mathbb{R}^{n\times r}$ whose columns span the orthogonal complement of the numerical kernel of $\widehat{Y}$. This $Q$ describes the minimal face of the PSD cone that contains $\widehat{Y}$ and supports the standard reparametrization $Y=QZQ^\mathsf{T}$ with $Z\succeq 0$. We give a high-level algorithm for the complete pipeline and then provide pseudocode for the main functions and helper subroutines, emphasizing how the code branches depending on the settings. We also explain how the “rational” route ($\mathrm{FF}=\mathbb{Q}$) differs from the algebraic number field route ($\mathrm{FF}=\mathbb{Q}(\alpha)$).



## Julia Codebase

The complete Julia codebase that we refer to as `rounding_modified.jl` is below.



````julia
using LinearAlgebra
using GenericLinearAlgebra
using Nemo
using RowEchelon


"""
    RoundingSettings(; kwargs...)

Settings for the single-matrix rounding pipeline.
Defaults are aligned with `src/rounding.jl`.
"""
struct RoundingSettings
    kernel_lll::Bool
    kernel_bits::Int64
    kernel_errbound::Float64
    kernel_round_errbound::Float64
    kernel_use_primal::Bool
    reduce_kernelvectors::Bool
    reduce_kernelvectors_cutoff::Int64
    reduce_kernelvectors_stepsize::Int64
    unimodular_transform::Bool
    approximation_decimals::Int64
    regularization::Float64
    normalize_transformation::Bool
    redundancyfactor::Int
    pseudo::Bool
    pseudo_columnfactor::Float64
    extracolumns_linindep::Bool
    function RoundingSettings(; kernel_lll=false,
                               kernel_bits=1000,
                               kernel_errbound=1e-10,
                               kernel_round_errbound=1e-15,
                               kernel_use_primal=true,
                               reduce_kernelvectors=true,
                               reduce_kernelvectors_cutoff=400,
                               reduce_kernelvectors_stepsize=200,
                               unimodular_transform=true,
                               approximation_decimals=40,
                               regularization=1e-20,
                               normalize_transformation=true,
                               redundancyfactor=10,
                               pseudo=true,
                               pseudo_columnfactor=1.05,
                               extracolumns_linindep=false)
        pseudo_columnfactor = max(pseudo_columnfactor, 1.0)
        new(kernel_lll,
            kernel_bits,
            kernel_errbound,
            kernel_round_errbound,
            kernel_use_primal,
            reduce_kernelvectors,
            reduce_kernelvectors_cutoff,
            reduce_kernelvectors_stepsize,
            unimodular_transform,
            approximation_decimals,
            regularization,
            normalize_transformation,
            redundancyfactor,
            pseudo,
            pseudo_columnfactor,
            extracolumns_linindep)
    end
end



"""
    _with_reduce_kernelvectors(settings, reduce_kernelvectors)

Return a copy of `settings` with only `reduce_kernelvectors` changed.
"""
function _with_reduce_kernelvectors(settings::RoundingSettings, reduce_kernelvectors::Bool)
    RoundingSettings(
        kernel_lll=settings.kernel_lll,
        kernel_bits=settings.kernel_bits,
        kernel_errbound=settings.kernel_errbound,
        kernel_round_errbound=settings.kernel_round_errbound,
        kernel_use_primal=settings.kernel_use_primal,
        reduce_kernelvectors=reduce_kernelvectors,
        reduce_kernelvectors_cutoff=settings.reduce_kernelvectors_cutoff,
        reduce_kernelvectors_stepsize=settings.reduce_kernelvectors_stepsize,
        unimodular_transform=settings.unimodular_transform,
        approximation_decimals=settings.approximation_decimals,
        regularization=settings.regularization,
        normalize_transformation=settings.normalize_transformation,
        redundancyfactor=settings.redundancyfactor,
        pseudo=settings.pseudo,
        pseudo_columnfactor=settings.pseudo_columnfactor,
        extracolumns_linindep=settings.extracolumns_linindep,
    )
end



"""
    _symmetrize_bigfloat(Yhat::AbstractMatrix)

Return a symmetric `BigFloat` copy of `Yhat`.

The result is `(Y + transpose(Y)) / 2`, where `Y = BigFloat.(Matrix(Yhat))`.
"""
function _symmetrize_bigfloat(Yhat::AbstractMatrix)
    size(Yhat, 1) == size(Yhat, 2) || error("Yhat must be square.")
    Y = BigFloat.(Matrix(Yhat))
    (Y + transpose(Y)) / 2
end

"""
    _identity_ff(FF, n::Int)

Return the `n × n` identity matrix over the field or ring `FF`.
"""
function _identity_ff(FF, n::Int)
    M = zero_matrix(FF, n, n)
    for i in 1:n
        M[i, i] = FF(1)
    end
    M
end

"""
    _to_real_matrix(Yhat::AbstractMatrix)

Return `(λmin, Y)` where `Y = _symmetrize_bigfloat(Yhat)` and `λmin` is the minimum
eigenvalue of `Symmetric(Y)`.
"""
function _to_real_matrix(Yhat::AbstractMatrix)
    Y = _symmetrize_bigfloat(Yhat)
    eigs = eigvals(Symmetric(Y))
    minimum(eigs), Y
end

"""
    generic_embedding(x, g; base_ring=BigFloat)

Embed exact elements into a real base ring using an approximation `g` of a field generator.
"""
function generic_embedding(x::AbsSimpleNumFieldElem, g; base_ring=BigFloat)
    F = parent(x)
    gb = base_ring(g)
    s = base_ring(0)
    for k = 0:degree(F)-1
        ck = coeff(x, k)
        s += generic_embedding(ck, g; base_ring=base_ring) * gb^k
    end
    s
end

"""
    generic_embedding(x::QQFieldElem, g; base_ring=BigFloat)

Embed a rational field element into `base_ring`.
"""
function generic_embedding(x::QQFieldElem, g; base_ring=BigFloat)
    base_ring(BigInt(numerator(x))) / base_ring(BigInt(denominator(x)))
end

"""
    generic_embedding(x::ZZRingElem, g; base_ring=BigFloat)

Embed an integer ring element into `base_ring`.
"""
function generic_embedding(x::ZZRingElem, g; base_ring=BigFloat)
    base_ring(BigInt(x))
end

"""
    generic_embedding(x::Rational{<:Integer}, g; base_ring=BigFloat)

Embed a Julia rational number into `base_ring`.
"""
function generic_embedding(x::Rational{<:Integer}, g; base_ring=BigFloat)
    base_ring(BigInt(numerator(x))) / base_ring(BigInt(denominator(x)))
end

@doc """
    generic_embedding(x::Integer, g; base_ring=BigFloat)

Embed an integer into `base_ring`.
"""
generic_embedding(x::Integer, g; base_ring=BigFloat) = base_ring(BigInt(x))

@doc """
    generic_embedding(x::AbstractFloat, g; base_ring=BigFloat)

Embed a floating-point value into `base_ring`.
"""
generic_embedding(x::AbstractFloat, g; base_ring=BigFloat) = base_ring(x)

"""
    _embed_vector(v, g)

Embed every entry of `v` into `BigFloat` using `generic_embedding` and generator
approximation `g`.
"""
function _embed_vector(v::AbstractVector, g)
    BigFloat[generic_embedding(x, g; base_ring=BigFloat) for x in v]
end


"""
    clindep(v, bits, errbound)

Find an integer relation `a` for vector `v` such that `sum(a[i]*v[i])` is below `errbound`.
"""
function clindep(v::AbstractVector, bits::Int, errbound)
    pow = 9
    F = AcbField(2^pow)
    u = F.(v)

    for p = 1:5:bits
        if p >= 2^(pow - 1)
            pow += 1
            F = AcbField(2^pow)
            u = F.(v)
        end
        a = ones(Int, length(u))
        try
            a = lindep(u, p)
        catch e
            if !(hasfield(typeof(e), :msg)) || !occursin("precision", e.msg)
                rethrow(e)
            end
        end
        err = abs(sum(u[i] * F(a[i]) for i in eachindex(u)))
        if err < errbound
            return a
        end
    end
    error("clindep failed to find a relation.")
end

"""
    clindep(v::AbstractMatrix, bits, errbound)

Find an integer relation among the rows of `v` with bounded residual.

This method returns integer coefficients `a` such that each entry of
`sum(a[i] * v[i, :])` is below `errbound` in absolute value.
"""
function clindep(v::AbstractMatrix, bits::Int, errbound)
    pow = 9
    F = AcbField(2^pow)
    u = map(F, transpose(v))
    for p = 1:5:bits
        if p >= 2^(pow - 1)
            pow += 1
            F = AcbField(2^pow)
            u = map(F, transpose(v))
        end
        a = ones(Int, size(u, 1))
        try
            a = lindep(u, p)
        catch e
            if !(hasfield(typeof(e), :msg)) || !occursin("precision", e.msg)
                rethrow(e)
            end
        end
        err = maximum(abs(sum(u[k, i] * F(a[i]) for i in eachindex(a))) for k in 1:size(u, 1))
        if err < errbound
            return a
        end
    end
    error("clindep failed to find a relation.")
end

"""
    roundx(x, g, deg; bits=100, errbound=1e-15)

Represent scalar `x` as coefficients in `QQ[g]/(minpoly)` up to degree `deg-1`.
Returns a vector of `deg` rational coefficients.
"""
function roundx(x::Number, g, deg::Int; bits=100, errbound=1e-15)
    start_vec = Any[x]
    for d = 0:deg-1
        push!(start_vec, g^d)
    end
    a = clindep(start_vec, bits, errbound)
    a[1] == 0 && error("Failed to recover scalar coefficient representation.")
    den = BigInt(a[1])
    [QQ(BigInt(-a[d+1]) // den) for d in 1:deg]
end

"""
    roundx(x::AbstractVector, g, deg; bits=100, errbound=1e-15)

Apply scalar `roundx` entrywise and return concatenated coefficient blocks.

If `deg = d` and `length(x) = n`, the output has length `d*n` and is ordered as
`[c₁(x₁), ..., c₁(xₙ), c₂(x₁), ..., c_d(xₙ)]`.
"""
function roundx(x::AbstractVector, g, deg::Int; bits=100, errbound=1e-15)
    coeff_blocks = [QQFieldElem[] for _ in 1:deg]
    for v in x
        coeffs = roundx(v, g, deg; bits=bits, errbound=errbound)
        for d in 1:deg
            push!(coeff_blocks[d], coeffs[d])
        end
    end
    vcat(coeff_blocks...)
end

"""
    _round_scalar_to_coeffs(x, g, deg; bits, errbound)

Recover degree-`deg` rational coefficients for scalar `x`, with fallback.

This first tries `roundx(x, g, deg; bits, errbound)`. If integer-relation recovery
fails, it falls back to rationalizing `x` into the constant coefficient.
"""
function _round_scalar_to_coeffs(x, g, deg::Int; bits::Int, errbound::Real)
    try
        return roundx(x, g, deg; bits=bits, errbound=errbound)
    catch
        q = rationalize(BigInt, BigFloat(x), tol=BigFloat(errbound))
        coeffs = fill(QQ(0), deg)
        coeffs[1] = QQ(q)
        return coeffs
    end
end


"""
    _select_independent_vectors(vectors, g; tol=BigFloat(1e-30))

Select a numerically linearly independent subsequence from `vectors`.

Independence is tested after embedding with `generic_embedding(·, g)` and applying
Gram-Schmidt with threshold `tol`.
"""
function _select_independent_vectors(vectors::Vector{<:AbstractVector}, g; tol=BigFloat(1e-30))
    chosen = eltype(vectors)[]
    orth = Vector{Vector{BigFloat}}()

    for v in vectors
        w = _embed_vector(v, g)
        for u in orth
            denom = dot(u, u)
            if denom != 0
                w .-= (dot(w, u) / denom) .* u
            end
        end
        if norm(w) > tol
            push!(chosen, v)
            push!(orth, w)
        end
    end
    chosen
end

"""
    _reduce_kernelvectors_qq(vectors)

Reduce rational kernel vectors by denominator clearing and LLL.

This function converts vectors to an integer lattice basis, applies `lll`, and
maps back to rational vectors.
"""
function _reduce_kernelvectors_qq(vectors::Vector{<:AbstractVector})
    isempty(vectors) && return vectors
    K = matrix(QQ, hcat(vectors...))
    for j in 1:ncols(K)
        l = lcm(Vector(denominator.(K[:, j])))
        K[:, j] *= l
    end
    Z = ZZ.(K)
    Zred = lll(transpose(Z))
    Kred = transpose(matrix(QQ, Zred))
    [Vector(Kred[:, j]) for j in 1:ncols(Kred)]
end

"""
    _complete_basis_with_standard_vectors(kernel_vectors, N, FF, g)

Complete `kernel_vectors` to a full basis of `FF^N` using standard vectors.

Candidates are accepted by numerical independence tests in the embedding defined by
`g`.
"""
function _complete_basis_with_standard_vectors(kernel_vectors, N::Int, FF, g)
    basis = copy(kernel_vectors)
    orth = Vector{Vector{BigFloat}}()

    for v in basis
        w = _embed_vector(v, g)
        for u in orth
            denom = dot(u, u)
            if denom != 0
                w .-= (dot(w, u) / denom) .* u
            end
        end
        if norm(w) > BigFloat(1e-30)
            push!(orth, w)
        end
    end

    for i in 1:N
        candidate_exact = fill(FF(0), N)
        candidate_exact[i] = FF(1)
        candidate_float = zeros(BigFloat, N)
        candidate_float[i] = 1

        for u in orth
            denom = dot(u, u)
            if denom != 0
                candidate_float .-= (dot(candidate_float, u) / denom) .* u
            end
        end

        if norm(candidate_float) > BigFloat(1e-30)
            push!(basis, candidate_exact)
            push!(orth, candidate_float)
        end
        length(basis) == N && break
    end
    basis
end

"""
    convert_system(FF, A, b)

Convert a linear system over a simple number field `FF` to a rational system.

Given `A * x = b` with entries in `FF = QQ(g)`, this expands coefficients in the
basis `(1, g, ..., g^(d-1))` and returns `(Atot, btot)` over `QQ`.
"""
function convert_system(FF, A, b)
    g = gen(FF)
    d = degree(FF)
    Ai = [matrix(QQ, coeff.(A, k)) for k=0:degree(FF)-1]
    btot = vcat([coeff.(b, k) for k in 0:d-1]...)

    n, m = size(A)
    Atot = zero_matrix(QQ, n * d, m * d)
    for i in 0:d-1
        for j in 0:d-1
            cur_gen = g^(i + j)
            for k in 0:d-1
                c = coeff(cur_gen, k)
                if c != 0
                    Atot[n * k + 1:n * (k + 1), m * j + 1:m * (j + 1)] += c * Ai[i + 1]
                end
            end
        end
    end
    return Atot, btot
end


"""
    detecteigenvectors(Yhat; FF=QQ, g=1, check_psd=true, settings=RoundingSettings())

Detect kernel vectors for one PSD matrix `Yhat` using high-precision SVD/RREF.

Returns vectors over field `FF`.
"""
function detecteigenvectors(Yhat::AbstractMatrix{T};
                            FF=QQ, g=1, check_psd=true,
                            settings::RoundingSettings=RoundingSettings()) where T
    min_eig, Y = _to_real_matrix(Yhat)
    if check_psd && min_eig < -10 * settings.kernel_errbound
        @warn "Yhat appears indefinite at tolerance $(settings.kernel_errbound). Minimum eigenvalue: $min_eig"
    end

    decomp = svd(Y)
    svals = decomp.S
    nzero = count(x -> abs(x) < settings.kernel_errbound, svals)
    if nzero == 0
        return Vector{Vector{typeof(FF(1))}}()
    end

    mat = Matrix(transpose(decomp.U[:, end-nzero+1:end]))
    RowEchelon.rref!(mat, settings.kernel_errbound)
    raw_vecs = [Vector(mat[i, :]) for i in axes(mat, 1) if norm(mat[i, :]) > settings.kernel_errbound]

    z = gen(FF)
    deg = degree(FF)
    kernel_vecs = Vector{typeof(FF(1))}[]
    for v in raw_vecs
        rounded = Vector{typeof(FF(1))}(undef, length(v))
        for i in eachindex(v)
            coeffs = _round_scalar_to_coeffs(v[i], g, deg; bits=settings.kernel_bits, errbound=settings.kernel_round_errbound)
            rounded[i] = sum(coeffs[k+1] * z^k for k = 0:deg-1)
        end
        res = Y * _embed_vector(rounded, g)
        if maximum(abs.(res)) <= 10 * settings.kernel_errbound
            push!(kernel_vecs, rounded)
        end
    end

    kernel_vecs
end


"""
    detecteigenvectors(Yhat, bits, errbound; FF=QQ, g=1)

Detect kernel vectors with the relation-matrix (`kernel_lll`) path.

This variant mirrors the original LLL-style approach: build field-power stacked
singular vectors, detect integer relations with `clindep`, recover a nullspace over
`ZZ`, and map back to vectors over `FF`. The start of the numerical kernel uses
the strict cutoff `|σ| < 1e-20` to match `rounding.jl`.
"""
function detecteigenvectors(Yhat::AbstractMatrix{T}, bits::Int, errbound=1e-15; FF=QQ, g=1) where T
    Y = _symmetrize_bigfloat(Yhat)
    tmp = svd(Y)
    svals = tmp.S
    M = tmp.U

    g_exact = gen(FF)
    deg = degree(FF)
    FtoQ = hcat([g^k * Matrix(I, size(M, 1), size(M, 1)) for k in 0:deg-1]...)
    FtoQ_exact = hcat([g_exact^k .* Matrix(I, size(M, 1), size(M, 1)) for k in 0:deg-1]...)
    Mstack = vcat([g^k .* M for k in 0:deg-1]...)

    list = Vector{ZZRingElem}[]
    first_zero = findfirst(x -> abs(x) < 1e-20, svals)
    if first_zero === nothing
        return Vector{Vector{typeof(FF(1))}}()
    end
    m = Mstack[:, first_zero:size(Mstack, 2)]

    if size(Y) == (1, 1) && abs(Y[1, 1]) <= 1e-6
        list = [ZZRingElem[ZZ(1)]]
    elseif size(m, 2) > 0
        A = zero_matrix(ZZ, 1, size(m, 1))
        s = collect(1:size(m, 1))
        while !isempty(s)
            l = clindep(m[s, :], bits, errbound)
            if deg == 1
                new_row = zero_matrix(ZZ, 1, size(m, 1))
                for (idx, val) in zip(s, l)
                    new_row[1, idx] = ZZ(val)
                end
                A = vcat(A, new_row)
            else
                lFF = reshape(FF.(l), :, 1)
                cur_A = FtoQ_exact[:, s] * lFF
                cur_b = [FF(0)]
                AQQ, _ = convert_system(FF, transpose(cur_A), cur_b)

                for r in eachrow(AQQ)
                    rowvals = collect(r)
                    lcm_row = lcm(denominator.(rowvals))
                    zrow = ZZ.(lcm_row .* rowvals)
                    new_row = zero_matrix(ZZ, 1, size(m, 1))
                    for (idx, val) in enumerate(zrow)
                        new_row[1, idx] = val
                    end
                    A = vcat(A, new_row)
                end
            end

            B = Nemo.rref(matrix(ZZ, A))[2]
            nonzero_rows = length([i for i in 1:nrows(B) if any(!iszero(B[i, j]) for j in 1:ncols(B))])
            if size(m, 1) - nonzero_rows - deg * size(m, 2) <= 0
                nmat = if size(B) == (1, 1) && B[1, 1] == 0
                    matrix(ZZ, 1, 1, [ZZ(1)])
                else
                    nullspace(B)[2]
                end
                list = [[nmat[k, i] for k in 1:nrows(nmat)] for i in 1:ncols(nmat)]
                break
            end

            all(iszero, l) && break
            idx = findfirst(x -> !iszero(x), l)
            idx === nothing && break
            deleteat!(s, idx)
        end
    end

    removals = Int[]
    for i in eachindex(list)
        res = Y * FtoQ * Vector{BigFloat}(list[i])
        if maximum(abs.(res)) > 1e-8
            push!(removals, i)
        end
    end
    deleteat!(list, removals)

    if length(list) != deg * size(m, 2)
        @warn "Not all kernel vectors detected by LLL path."
    end

    _select_independent_vectors([Vector(FtoQ_exact * v) for v in list], g; tol=BigFloat(errbound)^2)
end



"""
    simplify_kernelvectors(Yhat, finalvectors; FF=QQ, g=1, settings=RoundingSettings(), verbose=true)

Build a basis transformation `B = [V W]` from detected kernel vectors.

This follows the same reduction strategy as `src/rounding.jl`: RREF/nullspace/HNF
reduction for the default path and LLL-oriented reduction when `kernel_lll=true`.
"""
function simplify_kernelvectors(Yhat::AbstractMatrix, finalvectors;
                                FF=QQ, g=1,
                                settings::RoundingSettings=RoundingSettings(),
                                verbose=true)
    dm = _symmetrize_bigfloat(Yhat)
    isempty(finalvectors) && return _identity_ff(FF, size(dm, 1)), 0

    N = length(first(finalvectors))
    FF_kerneldim = length(finalvectors)
    if degree(FF) > 1
        lst = [vcat([coeff.(v .* gen(FF)^i, k) for k in 0:degree(FF)-1]...) for v in finalvectors for i in 0:degree(FF)-1]
    else
        lst = [QQ.(v) for v in finalvectors]
    end

    initial_maximum_number = max(
        maximum(maximum(denominator.(v)) for v in lst),
        maximum(maximum(abs.(numerator.(v))) for v in lst),
    )

    if settings.kernel_lll && settings.reduce_kernelvectors
        for i in eachindex(lst)
            lst[i] *= lcm(denominator.(lst[i]))
        end
        kernelvecs = transpose(matrix(ZZ, hcat(lst...)))
        B = transpose(lll(kernelvecs))
        kernel_dim = ncols(B)
        find_extra_vectors = true
    elseif settings.reduce_kernelvectors
        kernelvecs = transpose(matrix(QQ, hcat(lst...)))
        onehots = zeros(Int, nrows(kernelvecs))
        for i in 1:ncols(kernelvecs)
            success, j = all_except1(kernelvecs, i)
            if success && kernelvecs[j, i] == 1
                onehots[j] = i
            end
        end
        if any(==(0), onehots)
            error("The matrix of kernel vectors cannot easily be transformed to RREF.")
        end

        indices = unique(vcat(onehots, collect(1:ncols(kernelvecs))))
        indices_rev = [findfirst(==(k), indices) for k in 1:ncols(kernelvecs)]
        kernelvecs = kernelvecs[:, indices]

        if ncols(kernelvecs) > settings.reduce_kernelvectors_cutoff
            k = 1
            s = settings.reduce_kernelvectors_stepsize
            while true
                cols = unique(vcat(
                    collect(1:nrows(kernelvecs)),
                    collect(nrows(kernelvecs)+1:nrows(kernelvecs)+s*k),
                    collect(ncols(kernelvecs)-s*k+1:ncols(kernelvecs)),
                ))
                cols = [c for c in cols if 1 <= c <= ncols(kernelvecs)]
                part_kernelvecs = kernelvecs[:, cols]
                kernel_dim, B_part = reduction_step(part_kernelvecs; FF=FF, g=g)
                transform = transpose(B_part[1:kernel_dim, end-kernel_dim+1:end])
                reduced_kernelvecs = transform * kernelvecs

                if all(x -> denominator(x) == 1, reduced_kernelvecs)
                    reduced_kernelvecs = lll(ZZ.(reduced_kernelvecs))
                    B = transpose(reduced_kernelvecs)
                    find_extra_vectors = true
                    break
                else
                    for i in 1:nrows(reduced_kernelvecs)
                        l = lcm(denominator.(reduced_kernelvecs[i, :]))
                        reduced_kernelvecs[i, :] *= l
                    end
                    reduced_kernelvecs = lll(ZZ.(reduced_kernelvecs))
                    maxnum = maximum(abs.(reduced_kernelvecs))
                    verbose && @info "  The transform did not fully make the matrix integral. Maximum entry is $maxnum (compared to $initial_maximum_number)"
                    if maxnum <= initial_maximum_number
                        B = transpose(reduced_kernelvecs)
                        find_extra_vectors = true
                        break
                    else
                        k += 1
                    end
                end
            end
        else
            kernel_dim, B = reduction_step(kernelvecs; FF=FF, g=g)
            kernelvecs = transpose(B[:, end-kernel_dim+1:end])
            kernelvecs_reduced = transpose(lll(kernelvecs))
            if settings.unimodular_transform
                B[:, end-kernel_dim+1:end] = kernelvecs_reduced
                find_extra_vectors = false
            else
                B = kernelvecs_reduced
                find_extra_vectors = true
            end
        end
        B = B[indices_rev, :]
    else
        kernel_dim = degree(FF) * FF_kerneldim
        B = matrix(QQ, hcat(lst...))
        find_extra_vectors = true
    end

    final_maxnum = maximum(abs.(B))
    if degree(FF) > 1
        finalvectors = [sum(B[j * N + 1:(j + 1) * N, i] .* gen(FF)^j for j in 0:degree(FF)-1) for i in 1:ncols(B)]
        AF = ArbField(512)
        float_vecs = [generic_embedding.(x, AF(g), base_ring=AF) for x in finalvectors]

        if find_extra_vectors
            for i in 1:N
                vexact = zeros(FF, N)
                vexact[i] = FF(1)
                push!(finalvectors, vexact)
                v = zeros(AF, N)
                v[i] = AF(1)
                push!(float_vecs, v)
            end
        else
            finalvectors = [finalvectors[end-kernel_dim+1:end]..., finalvectors[1:end-kernel_dim]...]
            float_vecs = [float_vecs[end-kernel_dim+1:end]..., float_vecs[1:end-kernel_dim]...]
        end

        maxerror = maximum(maximum(abs.(dm * _embed_vector(v, g))) for v in finalvectors[1:kernel_dim])
        @assert maxerror < settings.kernel_errbound "The reduced kernel vectors are not kernel vectors (maximum error = $maxerror)"

        indices = Int[]
        orthogonalvectors = []
        for i in eachindex(float_vecs)
            for v in orthogonalvectors
                float_vecs[i] -= (dot(v, float_vecs[i]) // dot(v, v)) .* v
            end
            if dot(float_vecs[i], float_vecs[i]) > 1e-40
                push!(indices, i)
                push!(orthogonalvectors, float_vecs[i])
            end
            if length(orthogonalvectors) == N
                break
            end
        end
        finalvectors = finalvectors[indices]
        B = matrix(FF, hcat(finalvectors...))
    elseif find_extra_vectors
        maxerror = maximum(maximum(abs.(dm * BigFloat.(B[:, i]))) for i in 1:ncols(B))
        @assert maxerror < settings.kernel_errbound "The reduced kernel vectors are not kernel vectors (maximum error = $maxerror)"

        AF = ArbField(512)
        nlist = [Vector(AF.(B[:, i])) for i in 1:ncols(B)]
        for i in eachindex(nlist)
            for j in 1:i-1
                nlist[i] -= (dot(nlist[i], nlist[j]) // dot(nlist[j], nlist[j])) .* nlist[j]
            end
        end

        extravectors = zero_matrix(FF, N, N - ncols(B))
        j = 1
        for i in 1:N
            candidate = zeros(AF, N)
            candidate[i] = AF(1)
            for v in nlist
                candidate -= (dot(candidate, v) // dot(v, v)) .* v
            end
            if dot(candidate, candidate) > 1e-40
                push!(nlist, candidate)
                extravectors[i, j] = FF(1)
                j += 1
            end
            if length(nlist) == N
                break
            end
        end
        B = hcat(FF.(B), extravectors)
    else
        B = hcat(B[:, end-kernel_dim+1:end], B[:, 1:end-kernel_dim])
        maxerror = maximum(maximum(abs.(dm * BigFloat.(B[:, i]))) for i in 1:kernel_dim)
        @assert maxerror < settings.kernel_errbound "The reduced kernel vectors are not kernel vectors (maximum error = $maxerror)"
    end
    verbose && settings.reduce_kernelvectors && println("    After reduction, the maximum number in B is $final_maxnum")
    return B, FF_kerneldim
end



"""
    reduction_step(kernelvecs; FF=QQ, g=1)

Compute an integer basis transform from kernel rows in (near) RREF form.

Returns `(kernel_dim, B)`, where the last `kernel_dim` columns of `B` span the
kernel after denominator clearing and HNF-based transformation.
"""
function reduction_step(kernelvecs; FF=QQ, g=1)
    ns = transpose(nullspace_fromrref(kernelvecs)[2])
    for i in 1:nrows(ns)
        ns[i, :] *= lcm(Vector(denominator.(ns[i, :])))
    end
    kernel_dim, B = basis_nullspace_remaining(ZZ.(ns))
    return kernel_dim, B
end



"""
    basis_nullspace_remaining(x::ZZMatrix)

Split the transformed HNF basis into row-space and nullspace parts.

Returns `(k, B)` where `k` is the nullity and columns of `B` corresponding to the
nullspace define candidate kernel directions.
"""
function basis_nullspace_remaining(x::ZZMatrix)
    H, T = hnf_normalmultiplier_with_transform(transpose(x))
    for i in nrows(H):-1:1
        for j in 1:ncols(H)
            if !iszero(H[i, j])
                return nrows(H) - i, transpose(T)
            end
        end
    end
    return ncols(x), identity_matrix(x, ncols(x))
end



"""
    hnf_normalmultiplier_with_transform(A::ZZMatrix)

Return an HNF and a compatible left transformation multiplier.

For tall/square matrices this uses the augmented-matrix construction to recover a
stable multiplier that carries nullspace information.
"""
function hnf_normalmultiplier_with_transform(A::ZZMatrix)
    if nrows(A) < ncols(A)
        return hnf_with_transform(A)
    end
    mat = hcat(A, identity_matrix(ZZ, nrows(A)))
    H = hnf(mat)
    return H[:, 1:ncols(A)], H[:, ncols(A)+1:end]
end



"""
    nullspace_fromrref(M::QQMatrix)

Compute a rational nullspace basis, reusing `M` directly when it is already RREF.

Returns `(nullity, X)` where columns of `X` form a basis of `ker(M)`.
"""
function nullspace_fromrref(M::QQMatrix)
    m = nrows(M)
    n = ncols(M)
    if is_rref(M)
        A = M
        rank = 0
        for i in 1:m
            for j in 1:n
                if A[i, j] != 0
                    rank += 1
                    break
                end
            end
        end
    else
        rank, A = Nemo.rref(M)
    end

    nullity = n - rank
    R = base_ring(M)
    X = zero(M, n, nullity)
    if rank == 0
        for i in 1:nullity
            X[i, i] = one(R)
        end
    elseif nullity != 0
        pivots = zeros(Int, max(m, n))
        np = rank
        j = 1
        k = 1
        for i in 1:rank
            while is_zero_entry(A, i, j)
                pivots[np + k] = j
                j += 1
                k += 1
            end
            pivots[i] = j
            j += 1
        end
        while k <= nullity
            pivots[np + k] = j
            j += 1
            k += 1
        end
        for i in 1:nullity
            for j in 1:rank
                X[pivots[j], i] = -A[j, pivots[np + i]]
            end
            X[pivots[np + i], i] = one(R)
        end
    end
    return nullity, X
end

"""
    all_except1(m, k)

Check whether column `k` has exactly one nonzero entry.

Returns `(success, i)` where `success == true` means row `i` is the unique nonzero
position in that column.
"""
function all_except1(m, k)
    tot = 0
    indx = 0
    for i in 1:nrows(m)
        if m[i, k] != 0
            if tot == 1
                return false, indx
            end
            indx = i
            tot += 1
        end
    end
    if tot == 0
        return false, indx
    end
    return true, indx
end

"""
    basis_transformations(Yhat; FF=QQ, g=1, settings=RoundingSettings(), verbose=true)

Compute `(transpose(B), Binv, s)` for one PSD matrix `Yhat`, where `s` is the
detected nullity and columns `1:s` of `B` approximate a basis of `ker(Yhat)`.

If `settings.kernel_lll` is enabled, this first peels off near-zero diagonal
indices, runs the relation-based detector on the reduced submatrix, and lifts the
detected vectors back to ambient dimension before simplification.
"""
function basis_transformations(Yhat::AbstractMatrix;
                               FF=QQ, g=1,
                               settings::RoundingSettings=RoundingSettings(),
                               verbose=true)
    size(Yhat, 1) == size(Yhat, 2) || error("Yhat must be square.")
    Y = _symmetrize_bigfloat(Yhat)
    N = size(Y, 1)

    zerolist = Int[]
    list = Vector{typeof(FF(1))}[]
    if settings.kernel_lll
        for i in 1:N
            if abs(Y[i, i]) < settings.kernel_errbound
                push!(zerolist, i)
                v = fill(FF(0), N)
                v[i] = FF(1)
                push!(list, v)
            end
        end
    end
    nonzerolist = filter(x -> !(x in zerolist), collect(1:N))
    Ysub = Y[nonzerolist, nonzerolist]
    n = size(Ysub, 1)

    if n > 0
        if settings.kernel_lll
            prelist = detecteigenvectors(Ysub, settings.kernel_bits, settings.kernel_errbound; FF=FF, g=g)
            for l in prelist
                v = fill(FF(0), N)
                v[nonzerolist] = l
                push!(list, v)
            end
        else
            list = detecteigenvectors(Ysub; FF=FF, g=g, settings=settings)
        end

        if degree(FF) > 1
            coefficients = [coeff(x, k) for v in list for x in v for k in 0:degree(FF)-1]
        else
            coefficients = vcat(list...)
        end
        maxnum = maximum(abs.(numerator.(coefficients)), init=0)
        maxden = maximum(denominator.(coefficients), init=0)
        verbose && length(list) > 0 && println("    Matrix has $(length(list)) kernel vectors. Maximum numerator/denominator: $maxnum/$maxden")

        if settings.kernel_lll && degree(FF) > 1
            AF = ArbField(512)
            float_vecs = [generic_embedding.(v, AF(g), base_ring=AF) for v in list]
            orthogonalvectors = []
            lin_indep_list = Int[]
            for i in eachindex(float_vecs)
                for v in orthogonalvectors
                    float_vecs[i] -= (dot(v, float_vecs[i]) // dot(v, v)) .* v
                end
                if dot(float_vecs[i], float_vecs[i]) > 1e-40
                    push!(lin_indep_list, i)
                    push!(orthogonalvectors, float_vecs[i])
                end
            end
            verbose && println("  Finished GS to find linearly independent kernel vectors")
            list = list[lin_indep_list]
        end

        if length(list) > 0
            try
                B, num_kernelvecs = simplify_kernelvectors(Y, list; FF=FF, g=g, settings=settings, verbose=verbose)
            catch e
                if !settings.kernel_lll && settings.reduce_kernelvectors && e isa AssertionError
                    verbose && @warn "Kernel reduction failed in non-LLL path; retrying with reduce_kernelvectors=false."
                    fallback_settings = _with_reduce_kernelvectors(settings, false)
                    B, num_kernelvecs = simplify_kernelvectors(Y, list; FF=FF, g=g, settings=fallback_settings, verbose=verbose)
                else
                    rethrow(e)
                end
            end
            B = FF.(B)
        else
            num_kernelvecs = 0
            B = matrix(FF, Matrix(I, N, N))
        end
    else
        num_kernelvecs = length(list)
        if num_kernelvecs > 0
            B = matrix(FF, FF.(hcat(list...)))
        else
            B = matrix(FF, Matrix(I, N, N))
        end
    end
    Binv = Matrix(inv(B))

    if degree(FF) == 1 && settings.normalize_transformation
        lcms = [lcm(denominator.(Binv[i, :])) for i in 1:size(Binv, 1)]
        Binv = Diagonal(lcms) * Binv
        B = Matrix(B) * Diagonal([1 // x for x in lcms])
    else
        B = Matrix(B)
    end
    return transpose(B), Binv, num_kernelvecs
end

"""
    exact_solution(Yhat; FF=QQ, g=1, settings=RoundingSettings(), verbose=true)

Single-matrix version of `exact_solution`.
Returns `Q = B^{-T}[:, s+1:n]` and metadata, where `s = dim(ker(Yhat))`.

The returned rank fields are:
- `rank`: preserved compatibility alias (same as `rank_from_nullity`),
- `rank_from_nullity = n - s`,
- `rank_from_svd`: rank from high-precision singular values and `kernel_errbound`.
"""
function exact_solution(Yhat::AbstractMatrix;
                        FF=QQ, g=1,
                        settings::RoundingSettings=RoundingSettings(),
                        verbose=true)
    Bt, Binv, s = basis_transformations(Yhat; FF=FF, g=g, settings=settings, verbose=verbose)
    Q_exact = transpose(Binv)[:, s+1:end]
    Q = BigFloat.(generic_embedding.(Q_exact, g; base_ring=BigFloat))

    B = transpose(Bt)
    kernel_exact = s == 0 ? Matrix{eltype(B)}(undef, size(B, 1), 0) : B[:, 1:s]
    kernel_float = s == 0 ? zeros(BigFloat, size(B, 1), 0) : BigFloat.(generic_embedding.(kernel_exact, g; base_ring=BigFloat))
    Y = _symmetrize_bigfloat(Yhat)

    kernel_residual = s == 0 ? BigFloat(0) : maximum(abs.(Y * kernel_float))
    orthogonality = (size(Q, 2) == 0 || s == 0) ? BigFloat(0) : opnorm(transpose(kernel_float) * Q, Inf)
    rank_from_nullity = size(Q, 2)
    svals = svd(Y).S
    rank_from_svd = count(x -> abs(x) > settings.kernel_errbound, svals)
    rankY = rank_from_nullity

    return (
        Q = Q,
        Q_exact = Q_exact,
        rank = rankY,
        rank_from_nullity = rank_from_nullity,
        rank_from_svd = rank_from_svd,
        nullity = s,
        B = B,
        Bt = Bt,
        Binv = Binv,
        kernel_basis = kernel_float,
        diagnostics = (
            kernel_residual = kernel_residual,
            orthogonality = orthogonality,
            kernel_errbound = settings.kernel_errbound
        )
    )
end



"""
    example_single_matrix_pipeline()

Run the modified single-matrix pipeline on a small PSD matrix and return the result.
"""
function example_single_matrix_pipeline()
    A = [
        1.0 2.0 0.0 0.0;
        0.0 1.0 3.0 0.0
    ]
    Yhat = transpose(A) * A

    settings = RoundingSettings(
        kernel_errbound = 1e-12,
        kernel_round_errbound = 1e-14,
        kernel_lll = true,
        reduce_kernelvectors = true
    )

    result = exact_solution(Yhat; FF=QQ, g=1, settings=settings, verbose=true)
    return Yhat, result
end

if abspath(PROGRAM_FILE) == @__FILE__
    Yhat, result = example_single_matrix_pipeline()
    println("Input matrix size: ", size(Yhat))
    println("Detected rank: ", result.rank)
    println("Detected nullity: ", result.nullity)
    println("size(Q): ", size(result.Q))
    println("kernel residual: ", result.diagnostics.kernel_residual)
    println("orthogonality residual: ", result.diagnostics.orthogonality)
end
````



## Tutorial: Running `rounding_modified.jl` Step by Step

### 1.1 What we Discuss Here
This tutorial shows how to run the single-matrix rounding pipeline in
`rounding_modified.jl` step by step:

1. `detecteigenvectors`
2. `simplify_kernelvectors`
3. `basis_transformations`
4. `exact_solution`

We will see four examples:

1. Baseline step-by-step run on `Qx`.
2. Tolerance sweep on `Qx` (`1e-8` to `1e-10`).
3. Option comparison on `Qxnew` and `Zx` (`kernel_lll=false` vs `true`).
4. Number-field workflow (`FF = QQ(√2)`).

### 1.2 Setup and Imports
```julia
using Markdown
using LinearAlgebra
using Nemo

include("rounding_modified.jl")
```

### 1.3 Shared Helpers and Matrices
```julia

function inf_norm_identity_error(B, Binv, g=1)
    Bf = BigFloat.(generic_embedding.(B, g; base_ring=BigFloat))
    Binvf = BigFloat.(generic_embedding.(Binv, g; base_ring=BigFloat))
    n = size(B, 1)
    opnorm(Bf * Binvf - Matrix{BigFloat}(I, n, n), Inf)
end

Qx = [
     2 -1 0 0.217402307739 -0.332874917655 0 0.0010870115387 -0.00166437458827
    -1  1 0 -0.217402307739 0.332874917655 0.00192204561788 -0.0010870115387 0.00166437458827
     0  0 6.81292444806e-05 0.00317391904942 0.00273631052683 0 1.58695952471e-05 1.36815526341e-05
     0.217402307739 -0.217402307739 0.00317391904942 0.263591212134 0.0146877902575 0.00458359279312 0.00131795606067 7.34389512877e-05
    -0.332874917655 0.332874917655 0.00273631052683 0.0146877902575 0.291141352899 0.00307795438212 7.34389512877e-05 0.00145570676449
     0 0.00192204561788 0 0.00458359279312 0.00307795438212 2 2.29179639656e-05 1.53897719106e-05
     0.0010870115387 -0.0010870115387 1.58695952471e-05 0.00131795606067 7.34389512877e-05 2.29179639656e-05 6.58978030336e-06 3.67194756439e-07
    -0.00166437458827 0.00166437458827 1.36815526341e-05 7.34389512877e-05 0.00145570676449 1.53897719106e-05 3.67194756439e-07 7.27853382246e-06
]

Qxnew = [
     1.99999999983 -0.999999999451 0 0.693304542802 0.580759955609 0 0.00346652271401 0.00290379977804
    -0.999999999451 1.00000000037 0 -0.693304542802 -0.580759955609 0.00887012797498 -0.00346652271401 -0.00290379977804
     0 0 3.38157389103e-05 0.00305688283071 0.00293070879693 0 1.52844141536e-05 1.46535439846e-05
     0.693304542802 -0.693304542802 0.00305688283071 0.919482093467 0.623901099113 -0.00128065173054 0.00459741046272 0.00311950549556
     0.580759955609 -0.580759955609 0.00293070879693 0.623901099113 0.951647038028 -0.00387012797773 0.00311950549556 0.00475823518553
     0 0.00887012797498 0 -0.00128065173054 -0.00387012797773 1.99999999983 -6.40325865271e-06 -1.93506398886e-05
     0.00346652271401 -0.00346652271401 1.52844141536e-05 0.00459741046272 0.00311950549556 -6.40325865271e-06 2.2987975575e-05 1.55975274778e-05
     0.00290379977804 -0.00290379977804 1.46535439846e-05 0.00311950549556 0.00475823518553 -1.93506398886e-05 1.55975274778e-05 2.3792099189e-05
]

Zx = [
     1.72304247618e-05 -0.000485364354221 0.00103263537106 -0.000547943838735 -2.4268217711e-06 5.16317685528e-06 -2.73971919368e-06
    -0.000485364354221 2.93342907344 -3.43427915836 0.500723090255 0.0146671484315 -0.017171122879 0.00250402175276
     0.00103263537106 -3.43427915836 6.30764532962 -2.87305762642 -0.0171714760207 0.0315374460946 -0.0143659204743
    -0.000547943838735 0.500723090255 -2.87305762642 2.37223239156 0.00250366966739 -0.0143649038792 0.0118612682147
    -2.4268217711e-06 0.0146671484315 -0.0171714760207 0.00250366966739 7.33353267184e-05 -8.58560155395e-05 1.25203798444e-05
     5.16317685528e-06 -0.017171122879 0.0315374460946 -0.0143649038792 -8.58560155395e-05 0.000157683327705 -7.1827681107e-05
    -2.73971919368e-06 0.00250402175276 -0.0143659204743 0.0118612682147 1.25203798444e-05 -7.1827681107e-05 5.93068723582e-05
]



```

### 1.4 Example 1: Full Pipeline on `Qx` Step by Step
```julia
using Markdown
using LinearAlgebra
using Nemo

include("rounding_modified.jl")

Qx = [
     2 -1 0 0.217402307739 -0.332874917655 0 0.0010870115387 -0.00166437458827
    -1  1 0 -0.217402307739 0.332874917655 0.00192204561788 -0.0010870115387 0.00166437458827
     0  0 6.81292444806e-05 0.00317391904942 0.00273631052683 0 1.58695952471e-05 1.36815526341e-05
     0.217402307739 -0.217402307739 0.00317391904942 0.263591212134 0.0146877902575 0.00458359279312 0.00131795606067 7.34389512877e-05
    -0.332874917655 0.332874917655 0.00273631052683 0.0146877902575 0.291141352899 0.00307795438212 7.34389512877e-05 0.00145570676449
     0 0.00192204561788 0 0.00458359279312 0.00307795438212 2 2.29179639656e-05 1.53897719106e-05
     0.0010870115387 -0.0010870115387 1.58695952471e-05 0.00131795606067 7.34389512877e-05 2.29179639656e-05 6.58978030336e-06 3.67194756439e-07
    -0.00166437458827 0.00166437458827 1.36815526341e-05 7.34389512877e-05 0.00145570676449 1.53897719106e-05 3.67194756439e-07 7.27853382246e-06
]

settings = RoundingSettings(
    kernel_lll=false,
    kernel_errbound=1e-8,
    kernel_round_errbound=1e-12,
    reduce_kernelvectors=false
)

# Step 1: detecteigenvectors
# Use the non-LLL path first for robust nullity detection on this matrix.

vectors = detecteigenvectors(Qx; FF=QQ, g=1, settings=settings)
println("detected vectors: ", length(vectors))


# Step 2: simplify_kernelvectors
# Build `B = [V W]` where the first `s` columns represent kernel directions.

B, s = simplify_kernelvectors(Qx, vectors; FF=QQ, g=1, settings=settings, verbose=true)
println("kernel dimension s = ", s)
println("size(B) = ", size(B))


# Step 3: basis_transformations
# Compute `(Bt, Binv, s)` used by the final transformation.

Bt, Binv, s2 = basis_transformations(Qx; FF=QQ, g=1, settings=settings, verbose=true)
println("s2 = ", s2)
println("size(Bt) = ", size(Bt), ", size(Binv) = ", size(Binv))


# Step 4: exact_solution
# Extract `Q = B^{-T}[:, s+1:n]` and diagnostics.

out = exact_solution(Qx; FF=QQ, g=1, settings=settings, verbose=true)
println("size(Q) = ", size(out.Q))
println("nullity = ", out.nullity)
println("rank = ", out.rank)
println("rank_from_nullity = ", out.rank_from_nullity)
println("rank_from_svd = ", out.rank_from_svd)
println("kernel residual = ", out.diagnostics.kernel_residual)
println("orthogonality = ", out.diagnostics.orthogonality)

Bf = BigFloat.(generic_embedding.(out.B, 1; base_ring=BigFloat))
Binvf = BigFloat.(generic_embedding.(out.Binv, 1; base_ring=BigFloat))
println("||B*Binv - I||∞ = ", opnorm(Bf*Binvf - Matrix{BigFloat}(I, size(Qx,1), size(Qx,1)), Inf))
```

### 1.5 Example 2: Tolerance Sweep on `Qx` (`1e-8` to `1e-10`)
```julia
using Markdown
using Nemo

include("rounding_modified.jl")

Qx = [
     2 -1 0 0.217402307739 -0.332874917655 0 0.0010870115387 -0.00166437458827
    -1  1 0 -0.217402307739 0.332874917655 0.00192204561788 -0.0010870115387 0.00166437458827
     0  0 6.81292444806e-05 0.00317391904942 0.00273631052683 0 1.58695952471e-05 1.36815526341e-05
     0.217402307739 -0.217402307739 0.00317391904942 0.263591212134 0.0146877902575 0.00458359279312 0.00131795606067 7.34389512877e-05
    -0.332874917655 0.332874917655 0.00273631052683 0.0146877902575 0.291141352899 0.00307795438212 7.34389512877e-05 0.00145570676449
     0 0.00192204561788 0 0.00458359279312 0.00307795438212 2 2.29179639656e-05 1.53897719106e-05
     0.0010870115387 -0.0010870115387 1.58695952471e-05 0.00131795606067 7.34389512877e-05 2.29179639656e-05 6.58978030336e-06 3.67194756439e-07
    -0.00166437458827 0.00166437458827 1.36815526341e-05 7.34389512877e-05 0.00145570676449 1.53897719106e-05 3.67194756439e-07 7.27853382246e-06
]


# Tolerance sweep
# Smaller `kernel_errbound` is stricter. On near-singular matrices, this can change
detected nullity and final rank.

for eb in (1e-8, 1e-9, 1e-10)
    settings = RoundingSettings(
        kernel_lll=false,
        kernel_errbound=eb,
        kernel_round_errbound=1e-12,
        reduce_kernelvectors=false
    )
    vecs = detecteigenvectors(Qx; FF=QQ, g=1, settings=settings)
    out = exact_solution(Qx; FF=QQ, g=1, settings=settings, verbose=false)
    println("kernel_errbound = ", eb)
    println("  detecteigenvectors count = ", length(vecs))
    println("  nullity = ", out.nullity,
            ", rank_from_nullity = ", out.rank_from_nullity,
            ", rank_from_svd = ", out.rank_from_svd)
end
```

### 1.6 Example 3: `Qxnew` and `Zx` with `kernel_lll=false` vs `kernel_lll=true`
```julia
using Markdown
using Nemo

include("rounding_modified.jl")

Qxnew = [
     1.99999999983 -0.999999999451 0 0.693304542802 0.580759955609 0 0.00346652271401 0.00290379977804
    -0.999999999451 1.00000000037 0 -0.693304542802 -0.580759955609 0.00887012797498 -0.00346652271401 -0.00290379977804
     0 0 3.38157389103e-05 0.00305688283071 0.00293070879693 0 1.52844141536e-05 1.46535439846e-05
     0.693304542802 -0.693304542802 0.00305688283071 0.919482093467 0.623901099113 -0.00128065173054 0.00459741046272 0.00311950549556
     0.580759955609 -0.580759955609 0.00293070879693 0.623901099113 0.951647038028 -0.00387012797773 0.00311950549556 0.00475823518553
     0 0.00887012797498 0 -0.00128065173054 -0.00387012797773 1.99999999983 -6.40325865271e-06 -1.93506398886e-05
     0.00346652271401 -0.00346652271401 1.52844141536e-05 0.00459741046272 0.00311950549556 -6.40325865271e-06 2.2987975575e-05 1.55975274778e-05
     0.00290379977804 -0.00290379977804 1.46535439846e-05 0.00311950549556 0.00475823518553 -1.93506398886e-05 1.55975274778e-05 2.3792099189e-05
]

Zx = [
     1.72304247618e-05 -0.000485364354221 0.00103263537106 -0.000547943838735 -2.4268217711e-06 5.16317685528e-06 -2.73971919368e-06
    -0.000485364354221 2.93342907344 -3.43427915836 0.500723090255 0.0146671484315 -0.017171122879 0.00250402175276
     0.00103263537106 -3.43427915836 6.30764532962 -2.87305762642 -0.0171714760207 0.0315374460946 -0.0143659204743
    -0.000547943838735 0.500723090255 -2.87305762642 2.37223239156 0.00250366966739 -0.0143649038792 0.0118612682147
    -2.4268217711e-06 0.0146671484315 -0.0171714760207 0.00250366966739 7.33353267184e-05 -8.58560155395e-05 1.25203798444e-05
     5.16317685528e-06 -0.017171122879 0.0315374460946 -0.0143649038792 -8.58560155395e-05 0.000157683327705 -7.1827681107e-05
    -2.73971919368e-06 0.00250402175276 -0.0143659204743 0.0118612682147 1.25203798444e-05 -7.1827681107e-05 5.93068723582e-05
]


# Option comparison
# - `kernel_lll=false`: usually better for noisy near-kernel detection via `kernel_errbound`.
# - `kernel_lll=true`: uses strict internal onset (`|σ| < 1e-20`) for parity with `rounding.jl`,
#   so it may detect zero kernel vectors on near-singular numerical matrices.
# - The best `kernel_errbound` can depend on the matrix. Here we use matrix-specific
#   values chosen to keep the demo stable and illustrative.

for (name, Y, eb) in [("Qxnew", Qxnew, 1e-6), ("Zx", Zx, 5e-9)]
    println("\\nMatrix: ", name)
    println("  using kernel_errbound = ", eb)
    for (flag, label) in [(false, "non-LLL"), (true, "LLL")]
        settings = RoundingSettings(
            kernel_lll=flag,
            kernel_errbound=eb,
            kernel_round_errbound=1e-12,
            reduce_kernelvectors=false
        )
        try
            vecs = if flag
                detecteigenvectors(Y, settings.kernel_bits, settings.kernel_errbound; FF=QQ, g=1)
            else
                detecteigenvectors(Y; FF=QQ, g=1, settings=settings)
            end
            out = exact_solution(Y; FF=QQ, g=1, settings=settings, verbose=false)
            println("  ", label, ": detected=", length(vecs),
                    ", nullity=", out.nullity,
                    ", rank_from_svd=", out.rank_from_svd)
        catch err
            println("  ", label, ": run failed with ", typeof(err),
                    " (try adjusting kernel_errbound for this matrix)")
        end
    end
end
```

### 1.7 Example 4: Number Field Workflow (`FF = QQ(√2)`)
```julia
using Markdown
using Nemo

include("rounding_modified.jl")

R, x = polynomial_ring(QQ, "x")
K, a = number_field(x^2 - 2, "a")
g = sqrt(big(2))

Ysmall = [1.0 0.0; 0.0 0.0]

settings = RoundingSettings(
    kernel_lll=false,
    kernel_errbound=1e-10,
    kernel_round_errbound=1e-12,
    reduce_kernelvectors=false
)

md"""
### Number-field detector call
`g` is a numerical approximation of the generator of `K = QQ(√2)`.
"""
vectors = detecteigenvectors(Ysmall; FF=K, g=g, settings=settings)
println("detected vectors = ", length(vectors))
if !isempty(vectors)
    println("entry type = ", typeof(vectors[1][1]))
end


# Continue pipeline step-by-step in the number field

B, s = simplify_kernelvectors(Ysmall, vectors; FF=K, g=g, settings=settings, verbose=true)
println("s = ", s, ", size(B) = ", size(B))

Bt, Binv, s2 = basis_transformations(Ysmall; FF=K, g=g, settings=settings, verbose=true)
println("s2 = ", s2, ", size(Bt) = ", size(Bt), ", size(Binv) = ", size(Binv))

out = exact_solution(Ysmall; FF=K, g=g, settings=settings, verbose=true)
println("size(Q) = ", size(out.Q))
println("rank_from_nullity = ", out.rank_from_nullity, ", rank_from_svd = ", out.rank_from_svd)
```

### 1.8 Practical Guidance and Troubleshooting

#### Choosing options
1. Start with `kernel_lll=false`, `reduce_kernelvectors=false`, and a moderate
   `kernel_errbound` (for example `1e-8`).
2. If we need strict parity behavior with original LLL onset logic, use
   `kernel_lll=true`.
3. Enable `reduce_kernelvectors=true` when we want cleaner basis transforms and
   smaller coefficients, but keep an eye on diagnostics.

#### Interpreting diagnostics
1. `detected vectors = 0`:
   - In LLL mode, this can happen because of strict internal singular-value onset.
   - In non-LLL mode, try loosening `kernel_errbound`.
2. `rank_from_svd` vs `rank_from_nullity` mismatch:
   - This indicates tolerance/modeling sensitivity; revisit `kernel_errbound`.
3. Large kernel residual:
   - Tighten `kernel_round_errbound`.
   - Improve numerical quality of input matrix.

#### Recommended workflow for new matrices
1. Run Example 1 pattern with `kernel_lll=false`.
2. Sweep `kernel_errbound` as in Example 2.
3. Compare with `kernel_lll=true` as in Example 3.
4. Use number-field mode only when we explicitly want algebraic coefficients over a chosen field.




## Tutorial: Running `rounding_modified.jl` Step by Step

### 1.1 What we Discuss Here
This tutorial shows how to run the single-matrix rounding pipeline in
`rounding_modified.jl` step by step:

1. `detecteigenvectors`
2. `simplify_kernelvectors`
3. `basis_transformations`
4. `exact_solution`

We will see four examples:

1. Baseline step-by-step run on `Qx`.
2. Tolerance sweep on `Qx` (`1e-8` to `1e-10`).
3. Option comparison on `Qxnew` and `Zx` (`kernel_lll=false` vs `true`).
4. Number-field workflow (`FF = QQ(√2)`).

### 1.2 Setup and Imports
```julia
using Markdown
using LinearAlgebra
using Nemo

include("rounding_modified.jl")
```

### 1.3 Shared Helpers and Matrices
```julia

function inf_norm_identity_error(B, Binv, g=1)
    Bf = BigFloat.(generic_embedding.(B, g; base_ring=BigFloat))
    Binvf = BigFloat.(generic_embedding.(Binv, g; base_ring=BigFloat))
    n = size(B, 1)
    opnorm(Bf * Binvf - Matrix{BigFloat}(I, n, n), Inf)
end

Qx = [
     2 -1 0 0.217402307739 -0.332874917655 0 0.0010870115387 -0.00166437458827
    -1  1 0 -0.217402307739 0.332874917655 0.00192204561788 -0.0010870115387 0.00166437458827
     0  0 6.81292444806e-05 0.00317391904942 0.00273631052683 0 1.58695952471e-05 1.36815526341e-05
     0.217402307739 -0.217402307739 0.00317391904942 0.263591212134 0.0146877902575 0.00458359279312 0.00131795606067 7.34389512877e-05
    -0.332874917655 0.332874917655 0.00273631052683 0.0146877902575 0.291141352899 0.00307795438212 7.34389512877e-05 0.00145570676449
     0 0.00192204561788 0 0.00458359279312 0.00307795438212 2 2.29179639656e-05 1.53897719106e-05
     0.0010870115387 -0.0010870115387 1.58695952471e-05 0.00131795606067 7.34389512877e-05 2.29179639656e-05 6.58978030336e-06 3.67194756439e-07
    -0.00166437458827 0.00166437458827 1.36815526341e-05 7.34389512877e-05 0.00145570676449 1.53897719106e-05 3.67194756439e-07 7.27853382246e-06
]

Qxnew = [
     1.99999999983 -0.999999999451 0 0.693304542802 0.580759955609 0 0.00346652271401 0.00290379977804
    -0.999999999451 1.00000000037 0 -0.693304542802 -0.580759955609 0.00887012797498 -0.00346652271401 -0.00290379977804
     0 0 3.38157389103e-05 0.00305688283071 0.00293070879693 0 1.52844141536e-05 1.46535439846e-05
     0.693304542802 -0.693304542802 0.00305688283071 0.919482093467 0.623901099113 -0.00128065173054 0.00459741046272 0.00311950549556
     0.580759955609 -0.580759955609 0.00293070879693 0.623901099113 0.951647038028 -0.00387012797773 0.00311950549556 0.00475823518553
     0 0.00887012797498 0 -0.00128065173054 -0.00387012797773 1.99999999983 -6.40325865271e-06 -1.93506398886e-05
     0.00346652271401 -0.00346652271401 1.52844141536e-05 0.00459741046272 0.00311950549556 -6.40325865271e-06 2.2987975575e-05 1.55975274778e-05
     0.00290379977804 -0.00290379977804 1.46535439846e-05 0.00311950549556 0.00475823518553 -1.93506398886e-05 1.55975274778e-05 2.3792099189e-05
]

Zx = [
     1.72304247618e-05 -0.000485364354221 0.00103263537106 -0.000547943838735 -2.4268217711e-06 5.16317685528e-06 -2.73971919368e-06
    -0.000485364354221 2.93342907344 -3.43427915836 0.500723090255 0.0146671484315 -0.017171122879 0.00250402175276
     0.00103263537106 -3.43427915836 6.30764532962 -2.87305762642 -0.0171714760207 0.0315374460946 -0.0143659204743
    -0.000547943838735 0.500723090255 -2.87305762642 2.37223239156 0.00250366966739 -0.0143649038792 0.0118612682147
    -2.4268217711e-06 0.0146671484315 -0.0171714760207 0.00250366966739 7.33353267184e-05 -8.58560155395e-05 1.25203798444e-05
     5.16317685528e-06 -0.017171122879 0.0315374460946 -0.0143649038792 -8.58560155395e-05 0.000157683327705 -7.1827681107e-05
    -2.73971919368e-06 0.00250402175276 -0.0143659204743 0.0118612682147 1.25203798444e-05 -7.1827681107e-05 5.93068723582e-05
]



```

### 1.4 Example 1: Full Pipeline on `Qx` Step by Step
```julia
using Markdown
using LinearAlgebra
using Nemo

include("rounding_modified.jl")

Qx = [
     2 -1 0 0.217402307739 -0.332874917655 0 0.0010870115387 -0.00166437458827
    -1  1 0 -0.217402307739 0.332874917655 0.00192204561788 -0.0010870115387 0.00166437458827
     0  0 6.81292444806e-05 0.00317391904942 0.00273631052683 0 1.58695952471e-05 1.36815526341e-05
     0.217402307739 -0.217402307739 0.00317391904942 0.263591212134 0.0146877902575 0.00458359279312 0.00131795606067 7.34389512877e-05
    -0.332874917655 0.332874917655 0.00273631052683 0.0146877902575 0.291141352899 0.00307795438212 7.34389512877e-05 0.00145570676449
     0 0.00192204561788 0 0.00458359279312 0.00307795438212 2 2.29179639656e-05 1.53897719106e-05
     0.0010870115387 -0.0010870115387 1.58695952471e-05 0.00131795606067 7.34389512877e-05 2.29179639656e-05 6.58978030336e-06 3.67194756439e-07
    -0.00166437458827 0.00166437458827 1.36815526341e-05 7.34389512877e-05 0.00145570676449 1.53897719106e-05 3.67194756439e-07 7.27853382246e-06
]

settings = RoundingSettings(
    kernel_lll=false,
    kernel_errbound=1e-8,
    kernel_round_errbound=1e-12,
    reduce_kernelvectors=false
)

# Step 1: detecteigenvectors
# Use the non-LLL path first for robust nullity detection on this matrix.

vectors = detecteigenvectors(Qx; FF=QQ, g=1, settings=settings)
println("detected vectors: ", length(vectors))


# Step 2: simplify_kernelvectors
# Build `B = [V W]` where the first `s` columns represent kernel directions.

B, s = simplify_kernelvectors(Qx, vectors; FF=QQ, g=1, settings=settings, verbose=true)
println("kernel dimension s = ", s)
println("size(B) = ", size(B))


# Step 3: basis_transformations
# Compute `(Bt, Binv, s)` used by the final transformation.

Bt, Binv, s2 = basis_transformations(Qx; FF=QQ, g=1, settings=settings, verbose=true)
println("s2 = ", s2)
println("size(Bt) = ", size(Bt), ", size(Binv) = ", size(Binv))


# Step 4: exact_solution
# Extract `Q = B^{-T}[:, s+1:n]` and diagnostics.

out = exact_solution(Qx; FF=QQ, g=1, settings=settings, verbose=true)
println("size(Q) = ", size(out.Q))
println("nullity = ", out.nullity)
println("rank = ", out.rank)
println("rank_from_nullity = ", out.rank_from_nullity)
println("rank_from_svd = ", out.rank_from_svd)
println("kernel residual = ", out.diagnostics.kernel_residual)
println("orthogonality = ", out.diagnostics.orthogonality)

Bf = BigFloat.(generic_embedding.(out.B, 1; base_ring=BigFloat))
Binvf = BigFloat.(generic_embedding.(out.Binv, 1; base_ring=BigFloat))
println("||B*Binv - I||∞ = ", opnorm(Bf*Binvf - Matrix{BigFloat}(I, size(Qx,1), size(Qx,1)), Inf))
```

### 1.5 Example 2: Tolerance Sweep on `Qx` (`1e-8` to `1e-10`)
```julia
using Markdown
using Nemo

include("rounding_modified.jl")

Qx = [
     2 -1 0 0.217402307739 -0.332874917655 0 0.0010870115387 -0.00166437458827
    -1  1 0 -0.217402307739 0.332874917655 0.00192204561788 -0.0010870115387 0.00166437458827
     0  0 6.81292444806e-05 0.00317391904942 0.00273631052683 0 1.58695952471e-05 1.36815526341e-05
     0.217402307739 -0.217402307739 0.00317391904942 0.263591212134 0.0146877902575 0.00458359279312 0.00131795606067 7.34389512877e-05
    -0.332874917655 0.332874917655 0.00273631052683 0.0146877902575 0.291141352899 0.00307795438212 7.34389512877e-05 0.00145570676449
     0 0.00192204561788 0 0.00458359279312 0.00307795438212 2 2.29179639656e-05 1.53897719106e-05
     0.0010870115387 -0.0010870115387 1.58695952471e-05 0.00131795606067 7.34389512877e-05 2.29179639656e-05 6.58978030336e-06 3.67194756439e-07
    -0.00166437458827 0.00166437458827 1.36815526341e-05 7.34389512877e-05 0.00145570676449 1.53897719106e-05 3.67194756439e-07 7.27853382246e-06
]


# Tolerance sweep
# Smaller `kernel_errbound` is stricter. On near-singular matrices, this can change
detected nullity and final rank.

for eb in (1e-8, 1e-9, 1e-10)
    settings = RoundingSettings(
        kernel_lll=false,
        kernel_errbound=eb,
        kernel_round_errbound=1e-12,
        reduce_kernelvectors=false
    )
    vecs = detecteigenvectors(Qx; FF=QQ, g=1, settings=settings)
    out = exact_solution(Qx; FF=QQ, g=1, settings=settings, verbose=false)
    println("kernel_errbound = ", eb)
    println("  detecteigenvectors count = ", length(vecs))
    println("  nullity = ", out.nullity,
            ", rank_from_nullity = ", out.rank_from_nullity,
            ", rank_from_svd = ", out.rank_from_svd)
end
```

### 1.6 Example 3: `Qxnew` and `Zx` with `kernel_lll=false` vs `kernel_lll=true`
```julia
using Markdown
using Nemo

include("rounding_modified.jl")

Qxnew = [
     1.99999999983 -0.999999999451 0 0.693304542802 0.580759955609 0 0.00346652271401 0.00290379977804
    -0.999999999451 1.00000000037 0 -0.693304542802 -0.580759955609 0.00887012797498 -0.00346652271401 -0.00290379977804
     0 0 3.38157389103e-05 0.00305688283071 0.00293070879693 0 1.52844141536e-05 1.46535439846e-05
     0.693304542802 -0.693304542802 0.00305688283071 0.919482093467 0.623901099113 -0.00128065173054 0.00459741046272 0.00311950549556
     0.580759955609 -0.580759955609 0.00293070879693 0.623901099113 0.951647038028 -0.00387012797773 0.00311950549556 0.00475823518553
     0 0.00887012797498 0 -0.00128065173054 -0.00387012797773 1.99999999983 -6.40325865271e-06 -1.93506398886e-05
     0.00346652271401 -0.00346652271401 1.52844141536e-05 0.00459741046272 0.00311950549556 -6.40325865271e-06 2.2987975575e-05 1.55975274778e-05
     0.00290379977804 -0.00290379977804 1.46535439846e-05 0.00311950549556 0.00475823518553 -1.93506398886e-05 1.55975274778e-05 2.3792099189e-05
]

Zx = [
     1.72304247618e-05 -0.000485364354221 0.00103263537106 -0.000547943838735 -2.4268217711e-06 5.16317685528e-06 -2.73971919368e-06
    -0.000485364354221 2.93342907344 -3.43427915836 0.500723090255 0.0146671484315 -0.017171122879 0.00250402175276
     0.00103263537106 -3.43427915836 6.30764532962 -2.87305762642 -0.0171714760207 0.0315374460946 -0.0143659204743
    -0.000547943838735 0.500723090255 -2.87305762642 2.37223239156 0.00250366966739 -0.0143649038792 0.0118612682147
    -2.4268217711e-06 0.0146671484315 -0.0171714760207 0.00250366966739 7.33353267184e-05 -8.58560155395e-05 1.25203798444e-05
     5.16317685528e-06 -0.017171122879 0.0315374460946 -0.0143649038792 -8.58560155395e-05 0.000157683327705 -7.1827681107e-05
    -2.73971919368e-06 0.00250402175276 -0.0143659204743 0.0118612682147 1.25203798444e-05 -7.1827681107e-05 5.93068723582e-05
]


# Option comparison
# - `kernel_lll=false`: usually better for noisy near-kernel detection via `kernel_errbound`.
# - `kernel_lll=true`: uses strict internal onset (`|σ| < 1e-20`) for parity with `rounding.jl`,
#   so it may detect zero kernel vectors on near-singular numerical matrices.
# - The best `kernel_errbound` can depend on the matrix. Here we use matrix-specific
#   values chosen to keep the demo stable and illustrative.

for (name, Y, eb) in [("Qxnew", Qxnew, 1e-6), ("Zx", Zx, 5e-9)]
    println("\\nMatrix: ", name)
    println("  using kernel_errbound = ", eb)
    for (flag, label) in [(false, "non-LLL"), (true, "LLL")]
        settings = RoundingSettings(
            kernel_lll=flag,
            kernel_errbound=eb,
            kernel_round_errbound=1e-12,
            reduce_kernelvectors=false
        )
        try
            vecs = if flag
                detecteigenvectors(Y, settings.kernel_bits, settings.kernel_errbound; FF=QQ, g=1)
            else
                detecteigenvectors(Y; FF=QQ, g=1, settings=settings)
            end
            out = exact_solution(Y; FF=QQ, g=1, settings=settings, verbose=false)
            println("  ", label, ": detected=", length(vecs),
                    ", nullity=", out.nullity,
                    ", rank_from_svd=", out.rank_from_svd)
        catch err
            println("  ", label, ": run failed with ", typeof(err),
                    " (try adjusting kernel_errbound for this matrix)")
        end
    end
end
```

### 1.7 Example 4: Number Field Workflow (`FF = QQ(√2)`)
```julia
using Markdown
using Nemo

include("rounding_modified.jl")

R, x = polynomial_ring(QQ, "x")
K, a = number_field(x^2 - 2, "a")
g = sqrt(big(2))

Ysmall = [1.0 0.0; 0.0 0.0]

settings = RoundingSettings(
    kernel_lll=false,
    kernel_errbound=1e-10,
    kernel_round_errbound=1e-12,
    reduce_kernelvectors=false
)

md"""
### Number-field detector call
`g` is a numerical approximation of the generator of `K = QQ(√2)`.
"""
vectors = detecteigenvectors(Ysmall; FF=K, g=g, settings=settings)
println("detected vectors = ", length(vectors))
if !isempty(vectors)
    println("entry type = ", typeof(vectors[1][1]))
end


# Continue pipeline step-by-step in the number field

B, s = simplify_kernelvectors(Ysmall, vectors; FF=K, g=g, settings=settings, verbose=true)
println("s = ", s, ", size(B) = ", size(B))

Bt, Binv, s2 = basis_transformations(Ysmall; FF=K, g=g, settings=settings, verbose=true)
println("s2 = ", s2, ", size(Bt) = ", size(Bt), ", size(Binv) = ", size(Binv))

out = exact_solution(Ysmall; FF=K, g=g, settings=settings, verbose=true)
println("size(Q) = ", size(out.Q))
println("rank_from_nullity = ", out.rank_from_nullity, ", rank_from_svd = ", out.rank_from_svd)
```

### 1.8 Practical Guidance and Troubleshooting

#### Choosing options
1. Start with `kernel_lll=false`, `reduce_kernelvectors=false`, and a moderate
   `kernel_errbound` (for example `1e-8`).
2. If we need strict parity behavior with original LLL onset logic, use
   `kernel_lll=true`.
3. Enable `reduce_kernelvectors=true` when we want cleaner basis transforms and
   smaller coefficients, but keep an eye on diagnostics.

#### Interpreting diagnostics
1. `detected vectors = 0`:
   - In LLL mode, this can happen because of strict internal singular-value onset.
   - In non-LLL mode, try loosening `kernel_errbound`.
2. `rank_from_svd` vs `rank_from_nullity` mismatch:
   - This indicates tolerance/modeling sensitivity; revisit `kernel_errbound`.
3. Large kernel residual:
   - Tighten `kernel_round_errbound`.
   - Improve numerical quality of input matrix.

#### Recommended workflow for new matrices
1. Run Example 1 pattern with `kernel_lll=false`.
2. Sweep `kernel_errbound` as in Example 2.
3. Compare with `kernel_lll=true` as in Example 3.
4. Use number-field mode only when we explicitly want algebraic coefficients over a chosen field.

## Problem setting and goal

**Input.** A numerical symmetric matrix $\widehat{Y}\in\mathbb{S}^n$ that is expected to be PSD (or nearly PSD), typically returned by a floating-point SDP solver.

**Numerical kernel and rank (tolerance-based).** Fix a tolerance $\varepsilon>0$ (in the code: `settings.kernel_errbound`). We interpret $$K \;:=\; \operatorname{ker}(\widehat{Y}) \quad \text{numerically, as directions corresponding to singular values } \sigma < \varepsilon.$$ Let $s=\dim(K)$ be the numerical nullity and $r=n-s$ the corresponding numerical rank.

**Output (face basis).** The code returns a matrix $Q\in\mathbb{R}^{n\times r}$ whose columns form a basis of $K^\perp$. Equivalently, $\operatorname{range}(Q)=K^\perp \approx \operatorname{range}(\widehat{Y})$.

**Why $Q$ matters (minimal face / reparametrization).** If $\widehat{Y}\succeq 0$ and is rank-deficient, it lies on a proper face of the PSD cone. Knowing $K$ (or $K^\perp$) gives an explicit description of that face and enables dimension reduction and more robust exactification in larger rounding pipelines.

---

## Minimal faces of the PSD cone

Let $\mathbb{S}^n_{+}$ denote the cone of $n\times n$ real PSD matrices.

**Theorem (Minimal face of a PSD matrix).** Let $Y\in\mathbb{S}^n_{+}$ and let $K=\operatorname{ker}(Y)$ with $\dim(K)=s$ and $r=n-s$. Let $Q\in\mathbb{R}^{n\times r}$ have full column rank and satisfy $\operatorname{range}(Q)=K^\perp$. Then the minimal face of $\mathbb{S}^n_{+}$ containing $Y$ is $$\mathcal{F}(Y)=\{QZQ^\mathsf{T} : Z\in\mathbb{S}^r_+\}.$$ Equivalently, $$\mathcal{F}(Y)=\{X\succeq 0 : Xv=0 \ \forall v\in K\}.$$

**Interpretation.** The kernel $K$ imposes linear constraints $Xv=0$ for all $v\in K$. The parametrization $X=QZQ^\mathsf{T}$ automatically enforces these constraints and reduces dimension from $n\times n$ to $r\times r$.

**Connection to the code.** `rounding_modified.jl` focuses on constructing $Q$ from $\widehat{Y}$, by:
1.  finding exact (or exact-ish) kernel vectors spanning $K$ over an exact field $\mathrm{FF}$,  
2. constructing an invertible basis transformation $B=[V\ W]$ whose first columns $V$ span (approx.) $K$,  
3. outputting $Q = B^{-\mathsf{T}}[:,\,s+1:n]$, which spans $K^\perp$.

---

## The algebraic heart: why $Q = B^{-\mathsf{T}}_{(:,s+1:n)}$ spans $K^\perp$

The code constructs an invertible matrix $B\in\mathrm{FF}^{n\times n}$ and a number $s$ such that the first $s$ columns of $B$ form a basis of the kernel space $K$ (numerically).

Write the column partition $$B = [V\ W], \qquad V\in\mathrm{FF}^{n\times s},\;\; W\in\mathrm{FF}^{n\times (n-s)}.$$ Define $$Q_{\text{exact}} := B^{-\mathsf{T}}[:,\,s+1:n] \in \mathrm{FF}^{n\times (n-s)}.$$ After embedding $\mathrm{FF}\hookrightarrow \mathbb{R}$ numerically, the code returns $Q\in\mathbb{R}^{n\times (n-s)}$.

**Lemma (Orthogonality identity).** Let $B=[V\ W]\in\mathrm{FF}^{n\times n}$ be invertible. Then the blocks of $B^{-\mathsf{T}}=[V^\star\ W^\star]$ satisfy $V^\mathsf{T}W^\star = 0$. In particular, $\operatorname{range}(W^\star)\subseteq (\operatorname{range}(V))^\perp$ and $\dim(\operatorname{range}(W^\star))=n-s$.

**Proof.** From $B^\mathsf{T}B^{-\mathsf{T}}=I_n$ and the block form $B^\mathsf{T}=[V^\mathsf{T};\,W^\mathsf{T}]$, we get $$\begin{bmatrix} V^\mathsf{T}\\ W^\mathsf{T}\end{bmatrix}\begin{bmatrix}V^\star & W^\star\end{bmatrix}=\begin{bmatrix} V^\mathsf{T}V^\star & V^\mathsf{T}W^\star\\ W^\mathsf{T}V^\star & W^\mathsf{T}W^\star\end{bmatrix}=I_n.$$ Hence $V^\mathsf{T}W^\star = 0$. $\square$

**Conclusion.** If $V$ spans $K$ (the kernel of $\widehat{Y}$), then $W^\star=Q_{\text{exact}}$ spans $K^\perp$, which is exactly what we need to describe the minimal face.

---

## Exact arithmetic model: $\mathrm{FF}=\mathbb{Q}$ vs $\mathrm{FF}=\mathbb{Q}(\alpha)$

The pipeline manipulates vectors over an exact field $\mathrm{FF}$:
- **Rational route:** $\mathrm{FF}=\mathbb{Q}$ (Nemo’s `QQ`). Then $\deg(\mathrm{FF})=1$, and numbers are exact rationals.
- **Algebraic number field route:** $\mathrm{FF}=\mathbb{Q}(\alpha)$ where $\alpha$ is a root of a defining polynomial $p$. Then $\deg(\mathrm{FF})=d>1$. An element is represented as $\sum_{k=0}^{d-1} c_k z^k$ where $z=\mathrm{gen}(\mathrm{FF})$ and $c_k\in\mathbb{Q}$.

**Numerical embedding via an approximate generator.** To test residuals (e.g., whether a candidate vector is in $\operatorname{ker}(\widehat{Y})$), the code needs to embed exact elements into a real numeric ring (typically `BigFloat`): $$\phi_g:\ \mathrm{FF}\to\mathbb{R},\qquad \phi_g\!\left(\sum_{k=0}^{d-1} c_k z^k\right)=\sum_{k=0}^{d-1} c_k g^k,$$ where $g\in\mathbb{R}$ approximates the intended real root $\alpha$.

**Code mapping.** This is `generic_embedding(x,g;base_ring=BigFloat)` and `_embed_vector(v,g)`. For $\mathrm{FF}=\mathbb{Q}$, this reduces to standard rational-to-real conversion.

---

## Numerical-to-exact conversion via integer relations

A core operation is: given a floating scalar $x\in\mathbb{R}$ and a field generator approximation $g$, recover exact coefficients $(c_0,\dots,c_{d-1})\in\mathbb{Q}^d$ such that $$x \approx \sum_{k=0}^{d-1} c_k g^k.$$ This is an *integer relation* problem. The code uses Arb/Acb high-precision balls and `lindep`: find integers $a_0,\dots,a_d$ with $$a_0 x + a_1\cdot 1 + a_2 g + \cdots + a_d g^{d-1} \approx 0.$$ If $a_0\neq 0$, then $$x \approx -\sum_{k=1}^{d} \frac{a_k}{a_0} g^{k-1},\qquad\text{i.e., } c_{k-1} = -a_k/a_0\in\mathbb{Q}.$$

**Code mapping.**
- integer relation search: `clindep(...)`  
- coefficient recovery: `roundx(x,g,deg; bits, errbound)`  
- safe wrapper with fallback rationalization: `_round_scalar_to_coeffs`  

---

## Full pipeline: one-shot algorithm

We now present the pipeline at the level of `exact_solution(Yhat; FF,g,settings)`.

> **Algorithm (Single-matrix rounding / face detection: `exact_solution`).**
>
> - **Input:** symmetric numeric matrix $\widehat{Y}\in\mathbb{S}^n$ (expected PSD), field $\mathrm{FF}$, generator approximation $g\in\mathbb{R}$, settings $\mathsf{S}$  
> - **Output:** $Q\in\mathbb{R}^{n\times r}$ spanning $K^\perp$ where $K$ is the numerical kernel of $\widehat{Y}$  
>
> 1. $Y \gets (\widehat{Y}+\widehat{Y}^\mathsf{T})/2$ (symmetrize in high precision).  
> 2. $(B^\mathsf{T}, B^{-1}, s) \gets \texttt{BasisTransformations}(Y,\mathrm{FF},g,\mathsf{S})$ (with $s \approx \dim\operatorname{ker}(Y)$).  
> 3. $Q_{\text{exact}} \gets B^{-\mathsf{T}}[:,\,s+1:n]$ (exact basis for $K^\perp$ in $\mathrm{FF}$).  
> 4. $Q \gets \phi_g(Q_{\text{exact}})$ entrywise (embed to $\mathbb{R}$ via $g$).  
> 5. Compute diagnostics (kernel residual, orthogonality, ranks).  
> 6. Return $(Q, Q_{\text{exact}}, B, B^{-1}, s, \text{diagnostics})$.

**Code mapping.** This is implemented by: $$\texttt{exact\_solution} \;\to\; \texttt{basis\_transformations} \;\to\; \texttt{detecteigenvectors} \;\to\; \texttt{simplify\_kernelvectors}.$$

---

## `basis_transformations`: from kernel vectors to an invertible basis

The function `basis_transformations` performs three conceptual steps:
1. (optional) peel off obvious kernel coordinates by checking near-zero diagonal entries,  
2. detect kernel vectors by either the RREF path or the LLL path,  
3. reduce / complete these vectors into an invertible basis matrix $B=[V\ W]$ and return $(B^\mathsf{T},B^{-1},s)$.

> **Algorithm (`BasisTransformations`, maps to `basis_transformations`).**
>
> - **Input:** $Y\in\mathbb{S}^n$ (high precision symmetric), field $\mathrm{FF}$, generator approximation $g$, settings $\mathsf{S}$  
> - **Output:** $(B^\mathsf{T},B^{-1},s)$ with first $s$ columns of $B$ spanning numerical kernel  
>
> 1. Set $N\gets n$ and initialize $\mathcal{L}\gets[\ ]$ (list of exact kernel vectors in $\mathrm{FF}^N$).  
> 2. If $\mathsf{S}.\texttt{kernel\_lll}=\texttt{true}$:  
>    - $\mathcal{Z}\gets\{i:\left|Y_{ii}\right|<\mathsf{S}.\texttt{kernel\_errbound}\}$.  
>    - For each $i\in\mathcal{Z}$, add $e_i$ to $\mathcal{L}$ (standard basis vector, exact).  
>      Otherwise set $\mathcal{Z}\gets\emptyset$.  
> 3. $\mathcal{N}\gets \{1,\dots,N\}\setminus \mathcal{Z}$ and $Y_{\text{sub}} \gets Y[\mathcal{N},\mathcal{N}]$.  
> 4. If $\mathsf{S}.\texttt{kernel\_lll}=\texttt{true}$:  
>    - $\mathcal{K}_{\text{sub}} \gets \texttt{DetectEigenvectorsLLL}(Y_{\text{sub}},\mathsf{S}.\texttt{kernel\_bits},\mathsf{S}.\texttt{kernel\_errbound};\mathrm{FF},g)$.  
>    - For each $v\in\mathcal{K}_{\text{sub}}$, lift $v$ to ambient dimension by padding zeros on $\mathcal{Z}$ and placing entries on $\mathcal{N}$, then add lifted vector to $\mathcal{L}$.  
>      Else: $\mathcal{L}\gets \texttt{DetectEigenvectorsRREF}(Y_{\text{sub}};\mathrm{FF},g,\mathsf{S})$.  
> 5. If $\mathsf{S}.\texttt{kernel\_lll}=\texttt{true}$ and $\deg(\mathrm{FF})>1$, numerically filter $\mathcal{L}$ for linear independence (Gram–Schmidt in embedding).  
> 6. If $\mathcal{L}=\emptyset$, set $B\gets I_N$ and $s\gets 0$. Otherwise:  
>    - $(B,s)\gets \texttt{SimplifyKernelVectors}(Y,\mathcal{L};\mathrm{FF},g,\mathsf{S})$.  
>    - If non-LLL and the reduction assertion fails, retry with `reduce_kernelvectors=false` (implementation fallback).  
> 7. Compute $B^{-1}$ (matrix inverse over $\mathrm{FF}$).  
> 8. If $\deg(\mathrm{FF})=1$ and $\mathsf{S}.\texttt{normalize\_transformation}=\texttt{true}$: scale rows of $B^{-1}$ by row-wise LCMs of denominators and apply inverse scaling to $B$.  
> 9. Return $(B^\mathsf{T},B^{-1},s)$.

**Code mapping.**

- diagonal peeling: `zerolist` and adding standard vectors  
- kernel detection:  
  - `kernel_lll=false`: `detecteigenvectors(Ysub; ...)`  
  - `kernel_lll=true`: `detecteigenvectors(Ysub, bits, errbound; ...)`  
- simplification/reduction: `simplify_kernelvectors(Y,list; ...)`  
- normalization: the `lcms` diagonal scaling of `Binv` and reciprocal scaling of `B`  

---

## `detecteigenvectors`: two kernel detection modes

The pipeline supports two algorithmic paths for detecting kernel vectors.

### Mode A: SVD + RREF + entrywise rounding (default when `kernel_lll=false`)

Mathematically, the idea is:
1.  compute an orthonormal basis of the numerical left singular vectors corresponding to tiny singular values,  
2.  compute an RREF basis of that subspace (to obtain sparse/structured kernel vectors),  
3.  round each coordinate to the exact field $\mathrm{FF}$ via integer relations,  
4.  verify the residual $\left\|Yv\right\|$ is small (in the real embedding).

> **Algorithm (`DetectEigenvectorsRREF`, maps to `detecteigenvectors(Yhat; ...)`).**
>
> - **Input:** $Y\in\mathbb{S}^n$ (high precision), $\mathrm{FF}$, $g$, settings $\mathsf{S}$  
> - **Output:** list of exact kernel vectors $\mathcal{L}\subseteq \mathrm{FF}^n$  
>
> 1. Compute SVD: $Y = U\Sigma V^\mathsf{T}$.  
> 2. $s \gets \#\{i:\left|\sigma_i\right|<\mathsf{S}.\texttt{kernel\_errbound}\}$ (numerical nullity). If $s=0$, return empty list.  
> 3. $M \gets$ rows of $U^\mathsf{T}$ corresponding to the $s$ smallest singular values.  
> 4. Compute RREF of $M$ with tolerance $\mathsf{S}.\texttt{kernel\_errbound}$.  
> 5. Extract nonzero rows $\{w_1,\dots,w_s\}\subseteq \mathbb{R}^n$ (numerical kernel basis).  
> 6. Initialize $\mathcal{L}\gets [\ ]$. For each $j=1,\dots,s$:  
>    - For each coordinate $(w_j)_i$, recover exact coefficients in $\mathrm{FF}$ using integer relations.  
>    - Form $v_j \in \mathrm{FF}^n$.  
>    - If $\left\|Y\cdot \phi_g(v_j)\right\|_\infty \le C\cdot \mathsf{S}.\texttt{kernel\_errbound}$, add $v_j$ to $\mathcal{L}$.  
> 7. Return $\mathcal{L}$.

**Code mapping.**

- SVD: `decomp = svd(Y)` and `nzero = count(|σ| < kernel_errbound)`  
- RREF: `RowEchelon.rref!(mat, kernel_errbound)`  
- rounding: `_round_scalar_to_coeffs` + `rounded[i] = sum(coeffs[k+1]*z^k)`  
- residual check: `res = Y * _embed_vector(rounded,g)`  

**Implementation detail (important for users).** The current code *filters* rounded vectors using a residual threshold (`<= 10 * kernel_errbound`). This means the returned list may have fewer vectors than the SVD nullity estimate. (For a strict “fail-fast” behavior, one would error when a rounded candidate fails the residual test.)

### Mode B: LLL / integer relations on stacked singular vectors (when `kernel_lll=true`)

This mode is designed to be more robust when the numerical kernel space does not yield clean RREF vectors, especially in extension fields.

Idea (informal):
1. (a) compute singular vectors $u$ in the numerical kernel,  
2. (b) stack their multiples by powers of $g$ (or of the field generator),  
3. (c) find integer relations among these stacked rows using `lindep`,  
4. (d) build a relation matrix $A$ over $\mathbb{Z}$ and take its nullspace to recover kernel vectors.

For degree $d>1$, relations over $\mathrm{FF}$ are expanded into $\mathbb{Q}$-systems via coefficient expansion (`convert_system`).



> **Algorithm (`DetectEigenvectorsLLL`, maps to `detecteigenvectors(Yhat,bits,errbound;...)`).**
>
> - **Input:** $Y\in\mathbb{S}^n$ (high precision), $\mathrm{FF}$ of degree $d$, generator approx $g$, parameters `bits`, `errbound`  
> - **Output:** list of exact kernel vectors $\mathcal{L}\subseteq \mathrm{FF}^n$  
>
> 1. Compute SVD: $Y=U\Sigma V^\mathsf{T}$.  
> 2. Let $i_0$ be the first index with $\left|\sigma_{i_0}\right| < 10^{-20}$ (hard cutoff to match `rounding.jl`). If no such $i_0$, return empty list.  
> 3. $M \gets [U_{:,i_0},\dots,U_{:,n}]$ (numerical kernel singular vectors).  
> 4. Build stacked matrix $M_{\text{stack}} \gets \begin{bmatrix} M \\ gM \\ \vdots \\ g^{d-1}M \end{bmatrix}$.  
> 5. Initialize integer relation matrix $A\gets [0]$ over $\mathbb{Z}$.  
> 6. While not enough relations collected:  
>    - choose active row indices $s$ in $M_{\text{stack}}$  
>    - $l \gets \texttt{IntegerRelation}(M_{\text{stack}}[s,:], \texttt{bits}, \texttt{errbound})$ with $l\in\mathbb{Z}^{|s|}$  
>    - if $d=1$: append one row to $A$ encoding $l$ in ambient dimension  
>    - else: interpret $l$ as elements of $\mathrm{FF}$ and expand to $\mathbb{Q}$-constraints using coefficient expansion, then append expanded rows to $A$  
>    - reduce $A$ (e.g., RREF over $\mathbb{Z}$) and check if rank condition implies sufficient relations  
> 7. Compute integer nullspace basis $N$ of $A$ (columns of $N$ give integer combinations).  
> 8. Map each nullspace vector $v$ to $\mathrm{FF}^n$ via the stacked-basis map.  
> 9. Verify residuals and drop invalid candidates.  
> 10. Optionally select a numerically independent subset (via Gram–Schmidt embedding).  
> 11. Return kernel vectors in $\mathrm{FF}^n$.

**Code mapping.**
- stacked basis: `FtoQ`, `FtoQ_exact`, `Mstack`  
- relation finding: `clindep(m[s,:], bits, errbound)`  
- relation accumulation: `A = vcat(A,new_row)` and `B = Nemo.rref(matrix(ZZ,A))[2]`  
- nullspace extraction: `nullspace(B)[2]`  
- mapping back: `Vector(FtoQ_exact * v)`  

---

## `simplify_kernelvectors`: reduction, unimodularity, and basis completion

Given a list of exact kernel vectors $\mathcal{L}\subseteq \mathrm{FF}^n$, the function `simplify_kernelvectors` builds an invertible basis matrix $B\in\mathrm{FF}^{n\times n}$ whose first columns correspond to kernel directions (after reduction/reordering).

There are several subcases depending on options:
- `reduce_kernelvectors=false`: keep the provided kernel vectors (possibly large denominators) and complete to a full basis.  
- `reduce_kernelvectors=true`:  
  - if `kernel_lll=true`: clear denominators and apply LLL directly to shorten vectors;  
  - else: use an RREF-structured reduction (`reduction_step`) based on rational nullspaces and Hermite normal form (HNF), optionally keeping a unimodular transform.  
- `unimodular_transform`: if true, keep an ambient unimodular basis transform and replace only the kernel block by LLL-shortened vectors.

> **Algorithm (`SimplifyKernelVectors`, maps to `simplify_kernelvectors`).**
>
> - **Input:** $Y\in\mathbb{S}^n$, kernel vectors $\mathcal{L}\subseteq \mathrm{FF}^n$, settings $\mathsf{S}$  
> - **Output:** $(B,s)$ where $B\in\mathrm{FF}^{n\times n}$ invertible and first $s$ columns approximate $\operatorname{ker}(Y)$  
>
> 1. Convert $\mathcal{L}$ into a coefficient matrix over $\mathbb{Q}$ if $\deg(\mathrm{FF})>1$.  
> 2. Let $K$ be the matrix whose columns represent the kernel vectors (in $\mathbb{Q}$-coordinates).  
> 3. Set $s \gets \dim(\mathcal{L})$ (in $\mathrm{FF}$-vectors).  
> 4. If $\mathsf{S}.\texttt{kernel\_lll}=\texttt{true}$ and $\mathsf{S}.\texttt{reduce\_kernelvectors}=\texttt{true}$:  
>    - clear denominators so $K$ becomes integral  
>    - apply LLL to obtain shorter integral basis  
>    - set $B \gets$ the LLL-reduced basis (and later complete to full rank)  
> 5. Else if $\mathsf{S}.\texttt{reduce\_kernelvectors}=\texttt{true}$:  
>    - use RREF structure checks to find a stable pivot/column ordering (“onehot” heuristic)  
>    - if $ncols(K)$ is large: iteratively apply reduction on submatrices (controlled by cutoff/stepsize)  
>    - else: $(s_{\mathbb{Q}},B)\gets \texttt{ReductionStep}(K)$ (HNF-based transform), and optionally replace kernel block by LLL-shortened kernel rows while keeping unimodularity  
> 6. Else: set $B \gets$ matrix formed by the (possibly unreduced) kernel vectors; mark that extra vectors are needed.  
> 7. Complete $B$ to an invertible basis of $\mathrm{FF}^n$ (add standard vectors and select independent subset).  
> 8. Reorder columns so that kernel vectors occupy the first positions.  
> 9. Assert kernel residuals are within $\mathsf{S}.\texttt{kernel\_errbound}$.  
> 10. Return $(B,s)$.

**Code mapping.**
- coefficient expansion for $\deg(\mathrm{FF})>1$: `lst = [vcat(coeff.(...)) ...]`  
- direct LLL reduction in the (`kernel_lll` & `reduce`) case: `kernelvecs = transpose(matrix(ZZ,hcat(lst...)))` then `B = transpose(lll(kernelvecs))`  
- RREF/HNF reduction: `reduction_step` $\to$ `nullspace_fromrref` $\to$ `basis_nullspace_remaining` $\to$ `hnf_normalmultiplier_with_transform`  
- completing the basis: the `find_extra_vectors` logic and Gram–Schmidt independence screening  

---

## Key helper subroutines (pseudocode + code mapping)

This section summarizes the most important helper routines and their mathematical meaning.

### Nullspace from RREF over $\mathbb{Q}$ (`nullspace_fromrref`)

Given an RREF matrix $M\in\mathbb{Q}^{m\times n}$, the nullspace basis can be written explicitly by choosing free variables. The code uses a standard pivot/free-variable construction.

> **Algorithm (`NullspaceFromRREF`, maps to `nullspace_fromrref`).**
>
> - **Input:** $M\in\mathbb{Q}^{m\times n}$, optionally already in RREF  
> - **Output:** $(\operatorname{null}(M),X)$ where columns of $X\in\mathbb{Q}^{n\times \operatorname{null}(M)}$ span $\operatorname{ker}(M)$  
>
> 1. If $M$ is not in RREF, compute $(r,A)\gets \texttt{RREF}(M)$; else set $A\gets M$ and $r\gets \operatorname{rank}(M)$.  
> 2. $\ell \gets n-r$ (nullity).  
> 3. Identify pivot columns $p_1,\dots,p_r$ and free columns $f_1,\dots,f_\ell$.  
> 4. Initialize $X\gets 0\in\mathbb{Q}^{n\times \ell}$.  
> 5. For $j=1,\dots,\ell$:  
>    - set $X[f_j,j]=1$  
>    - for $i=1,\dots,r$: set $X[p_i,j] = -A[i,f_j]$  
> 6. Return $(\ell,X)$.

**Code mapping.** This is almost line-for-line the implementation in `nullspace_fromrref(M::QQMatrix)`.

### HNF-based transformation that exposes nullspace (`hnf_normalmultiplier_with_transform`)

Hermite normal form (HNF) computes an integer row reduction. The augmented construction $$\left[\;A\;\middle|\;I\;\right] \xrightarrow{\text{HNF}} \left[\;H\;\middle|\;T\;\right]$$ produces a transformation matrix $T$ that encodes (among other things) nullspace structure.

> **Algorithm (`HNFNormalMultiplierWithTransform`, maps to `hnf_normalmultiplier_with_transform`).**
>
> - **Input:** integer matrix $A\in\mathbb{Z}^{m\times n}$  
> - **Output:** $(H,T)$ such that $T$ is a left multiplier and $H$ is an HNF-like reduction  
>
> 1. If $m<n$: compute $(H,T)\gets \texttt{HNFWithTransform}(A)$ (unique in this regime).  
> 2. Else: form augmented matrix $[A\;|\;I_m]$, compute HNF $[A\;|\;I_m]\mapsto [H\;|\;T]$.  
> 3. Return $(H,T)$.

**Code mapping.** The augmented strategy is implemented by `mat = hcat(A, identity_matrix(ZZ,nrows(A)))`, then `H = hnf(mat)`, then splitting `H[:,1:ncols(A)]` and `H[:,ncols(A)+1:end]`.

### Reduction step (`reduction_step` and `basis_nullspace_remaining`)

This step turns a row-space description (kernel vectors in near-RREF form) into a basis transformation. The conceptual steps are:
1. (i) compute a rational nullspace basis of the kernel-vector matrix,  
2. (ii) clear denominators to obtain an integer matrix,  
3. (iii) use HNF to get a transformation whose last columns represent the kernel directions nicely.

> **Algorithm (`ReductionStep`, maps to `reduction_step` + `basis_nullspace_remaining`).**
>
> - **Input:** kernel matrix $K\in\mathbb{Q}^{s\times n}$ (rows represent kernel constraints)  
> - **Output:** $(s_{\mathbb{Z}},B)$ such that the last $s_{\mathbb{Z}}$ columns of $B$ encode kernel directions  
>
> 1. $(\ell,X)\gets \texttt{NullspaceFromRREF}(K)$.  
> 2. $N\gets X^\mathsf{T}$ (rows will be cleared of denominators).  
> 3. For each row $i$ of $N$, multiply row $i$ by $\mathrm{lcm}$ of denominators to make it integral.  
> 4. $(H,T)\gets \texttt{HNFNormalMultiplierWithTransform}(N^\mathsf{T})$.  
> 5. Determine how many trailing columns of $T^\mathsf{T}$ correspond to nullspace directions.  
> 6. Return $(s_{\mathbb{Z}},\,B=T^\mathsf{T})$.

**Code mapping.** See: `ns = transpose(nullspace_fromrref(kernelvecs)[2])`, denominator clearing, then `basis_nullspace_remaining(ZZ.(ns))`.

---

## Option-dependent paths: what changes when settings change?

The pipeline behavior is controlled primarily by $$\texttt{kernel\_lll},\quad \texttt{kernel\_bits},\quad \texttt{kernel\_errbound},\quad \texttt{kernel\_round\_errbound},\quad \texttt{reduce\_kernelvectors},\quad \texttt{unimodular\_transform},\quad \texttt{normalize\_transformation}.$$

### Kernel detection switches

- **`kernel_lll=false`:** use Mode A (SVD $\to$ RREF $\to$ entrywise rounding). Typically faster and simpler when the numerical kernel is clean.  
- **`kernel_lll=true`:** use Mode B (integer relations / LLL-style). More robust in tricky cases, and closer to the original `rounding.jl` “relation matrix” approach.

### Kernel reduction switches

- **`reduce_kernelvectors=false`:** do not try to shorten vectors; simply complete them to a basis. This can yield large denominators and large entries in $B$.  
- **`reduce_kernelvectors=true`:** apply LLL/HNF-based reduction heuristics to keep numbers smaller and make $B$ better conditioned.

### Unimodularity and basis transformation quality

- **`unimodular_transform=true`** (in the non-LLL reduction path): keep a unimodular ambient transform and only replace the kernel block with LLL-shortened vectors. This tends to preserve invertibility and keeps transformations “integral” in a strong sense.  
- **`unimodular_transform=false`:** use the shortened kernel vectors directly and then complete the basis; unimodularity is not preserved.

### Normalization (rational-only)

When $\mathrm{FF}=\mathbb{Q}$ and `normalize_transformation=true`, the code rescales the inverse matrix $B^{-1}$ by row-wise least common multiples of denominators so that subsequent transformed objects would avoid unnecessary rational denominators. This does not change the span decomposition; it only rescales the basis.

---

## Rational route vs algebraic number field route: conceptual differences

### What is being recovered?

- **Rational route ($\mathrm{FF}=\mathbb{Q}$):** each rounded entry is a rational number $p/q$. Kernel vectors lie in $\mathbb{Q}^n$.  
- **Number field route ($\mathrm{FF}=\mathbb{Q}(\alpha)$, $\deg=d>1$):** each rounded entry is an element $\sum_{k=0}^{d-1} c_k z^k$. Kernel vectors lie in $\mathrm{FF}^n$.

### How rounding works (entrywise)

- **$\mathrm{FF}=\mathbb{Q}$:** effectively recover $x\approx c_0$ (degree $1$). Integer relation search reduces to a rational reconstruction problem.  
- **$\mathrm{FF}=\mathbb{Q}(\alpha)$:** recover $x\approx \sum_{k=0}^{d-1} c_k g^k$ where $g\approx\alpha$. This requires integer relation detection in dimension $d+1$ per scalar entry (one unknown relation coefficient per basis element plus $x$).

### How reduction differs

- In $\mathbb{Q}$, reduction is performed directly with rational/integer linear algebra (LLL, HNF) on coefficient matrices.  
- In $\mathbb{Q}(\alpha)$, vectors are first expanded into $\mathbb{Q}$-coefficients (size blow-up by factor $d$), then reduced, then lifted back into $\mathrm{FF}$.

### Why a real embedding is essential in the number field case

To test whether a purported exact vector $v\in\mathrm{FF}^n$ is truly in the kernel of $\widehat{Y}$, the code evaluates $v$ numerically via $g$ and checks $\widehat{Y}\cdot \phi_g(v)\approx 0$. This is implemented by `generic_embedding` and `_embed_vector`.

---

## Diagnostics returned by `exact_solution`

The pipeline returns several useful metrics:
- **kernel residual:** $\max\left|Y\,V\right|$ where $V$ is the returned kernel basis (embedded numerically). This checks whether the detected kernel directions are indeed near-kernel.  
- **orthogonality residual:** $\left\|V^\mathsf{T}Q\right\|_\infty$. This checks that $Q$ is (numerically) orthogonal to the kernel basis.  
- **rank estimates:**  
  - `rank_from_nullity`: $n-s$,  
  - `rank_from_svd`: number of singular values $\sigma > \texttt{kernel\_errbound}$.

---

## Practical guidance on choosing settings (brief)

- Start with `kernel_lll=false` for speed. If kernel detection is unstable, try `kernel_lll=true`.  
- Tighten `kernel_errbound` if you observe spurious kernel vectors; loosen it if true kernel vectors are missed.  
- For exactification difficulty (especially in number fields), increase `kernel_bits` and/or relax `kernel_round_errbound`.  
- If $B$ becomes huge or ill-conditioned, keep `reduce_kernelvectors=true` and consider `unimodular_transform=true`.

---

## Summary

The single-matrix pipeline in `rounding_modified.jl` implements a mathematically standard face-detection routine for PSD matrices:
1. detect the numerical kernel of $\widehat{Y}$ and recover exact kernel vectors in an exact field $\mathrm{FF}$,  
2. reduce/complete these vectors into an invertible basis matrix $B=[V\ W]$ with $V$ spanning the kernel,  
3. output $Q=B^{-\mathsf{T}}_{(:,s+1:n)}$ which spans $K^\perp$ and describes the minimal face $\{QZQ^\mathsf{T}\}$.

The rational route ($\mathrm{FF}=\mathbb{Q}$) yields exact rational bases, while the algebraic number field route ($\mathrm{FF}=\mathbb{Q}(\alpha)$) uses integer-relation methods to recover algebraic coefficients and a real embedding to validate and post-process the results.

