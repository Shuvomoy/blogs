@def title = "How to solve semidefinite optimization problems in Julia"
@def published ="January 1, 2021"
@def tags =["programming", "Julia", "optimization"]

# How to solve semidefinite optimization problems in Julia
**Shuvomoy Das Gupta**

*January 1, 2021*

In this blog, we discuss how to solve semidefinite programs (SDPs) in ``Julia`` using ``Convex.jl``. We will consider two examples: (i) standard form semidefinite program, (ii) a somewhat more complicated semidefinite program.

---

**Table of contents**

\toc

---

### Standard form sdp

We consider optimization problem of the form:

$$
\begin{align*}
\begin{array}{ll}
\textrm{minimize} & \mathbf{trace}(CX)\\
\textrm{subject to} & \mathbf{trace}(A_{i}X)=b_{i},\\
 & X\succeq0,
\end{array} & i=1,\ldots,m
\end{align*}
$$


where $X\in\mathbf{S}^{n}$ is the decision variable, and each of the $A_{i}$ matrices and $C$ are also in $\mathbf{S}^{n}$. By the notation $\mathbf{S}^{n}$, we denote the set of all symmetric $n\times n$ matrices.

First, we load the necessary `Julia` packages.


```julia
using SCS, COSMO, MosekTools, JuMP, LinearAlgebra

using BenchmarkTools
```

Let us create data $A,C,b$ randomly.


```julia
function random_mat_create(n)
    # this function creates a symmetric n×n matrix
    A = randn(n,n)
    A = A'*A
    A = (A+A')/2
    return A
end
```

Here is the data generation process, please change it to your need.


```julia
n = 10
m = 20
# set of all data matrices A_i
# the data matrix A = [A1 A2 A3 ....]
A = zeros(n, m*n)
b = zeros(m)
# just ensuring our problem is feasible
X_test = rand(n,n)
X_test = X_test'*X_test
X_test = (X_test+X_test')/2
for i in 1:m
    A[:, (i-1)*n+1:i*n] .= random_mat_create(n)
    b[i] = tr(A[:, (i-1)*n+1:i*n]*X_test)
end
C = abs.(random_mat_create(n))
```

The following function solves the underlying SDP.


```julia

function solve_SDP(A, b, C; solver_name=:COSMO)

# Create variable
    if solver_name == :COSMO
        model = Model(COSMO.Optimizer)
    elseif solver_name == :Mosek
        model = Model(optimizer_with_attributes(Mosek.Optimizer))
    end

    set_silent(model)

    @variable(model, X[1:n, 1:n], PSD)


    @objective(model, Min, tr(C * X));
    for j in 1:m
        A_j = A[:, (j - 1) * n + 1:j * n]
        @constraint(model, tr(A_j * X) == b[j])
    end

    optimize!(model)

    status = JuMP.termination_status(model)
    X_sol = JuMP.value.(X)
    obj_value = JuMP.objective_value(model)

    return status, X_sol, obj_value

end
```

Time to solve the problem.


```julia
status, X_sol, obj_value = solve_SDP(A, b, C; solver_name=:Mosek)
```


```julia
out: (MathOptInterface.OPTIMAL, [2.907044311952373 1.7130367276575142 … -0.056145513617222656 3.0230926674218024; 1.7130367276575142 3.419039624378557 … 1.0871948703965775 2.0577919984154334; … ; -0.056145513617222656 1.0871948703965775 … 0.7127834057861842 0.5195987956934747; 3.0230926674218024 2.0577919984154334 … 0.5195987956934747 4.234108180247669], 939.2696385581793)
```

Lets see which solver is faster, `COSMO` or `Mosek`.


```julia
b1 = @benchmark solve_SDP(A, b, C; solver_name=:COSMO)

println("benchmark for COSMO")
println("*************************")
io = IOBuffer()
show(io, "text/plain", b1)
s = String(take!(io))
println(s)
```


```julia
benchmark for COSMO
*************************
BenchmarkTools.Trial:
  memory estimate:  4.33 MiB
  allocs estimate:  35786
  --------------
  minimum time:     36.964 ms (0.00% GC)
  median time:      40.552 ms (0.00% GC)
  mean time:        40.787 ms (1.14% GC)
  maximum time:     46.920 ms (9.97% GC)
  --------------
  samples:          123
  evals/sample:     1
```


```julia
b2 = @benchmark solve_SDP(A, b, C; solver_name=:Mosek)

println("benchmark for Mosek")
println("***************************")
io = IOBuffer()
show(io, "text/plain", b2)
s = String(take!(io))
println(s)
```


```julia
benchmark for Mosek
***************************
BenchmarkTools.Trial:
  memory estimate:  3.81 MiB
  allocs estimate:  32015
  --------------
  minimum time:     6.647 ms (0.00% GC)
  median time:      7.524 ms (0.00% GC)
  mean time:        8.374 ms (5.00% GC)
  maximum time:     19.410 ms (25.35% GC)
  --------------
  samples:          597
  evals/sample:     1
```

So, on average, `Mosek` seems to be 5 times faster than `COSMO`.

### Complicated sdp

Denote,
$$
\begin{align*}
J & =\{(i,j)\mid j=i+1,i\in[0:N-1]\}\cup\{(i,j)\mid i=\star,j\in[0:N]\},
\end{align*}
$$
 and
$$
\mathbf{w}_{0}=e_{1}\in\mathbf{R}^{N+2},\mathbf{g}_{i}=e_{i+2}\in\mathbf{R}^{N+2},\mathbf{f}_{i}=e_{i+1}\in\mathbf{R}^{N+1},
$$
where $e_{i}$ is the unit vector with $i$th component equal to $1$ and the rest being zero. Next, define
$$
\begin{align*}
 & S\left(\tau,\{\lambda_{i,j}\},\{\alpha_{i,j}^{\prime}\}\right)\\
= & c_{w}\tau\mathbf{w}_{0}\mathbf{w}_{0}^{\top}+\\
 & \frac{1}{2L}\Bigg[\sum_{i\in[0:N-1]}\lambda_{i,i+1}\left\{ (\mathbf{g}_{i}-\mathbf{g}_{i+1})\odot(\mathbf{g}_{i}-\mathbf{g}_{i+1})\right\} +\\
 & \sum_{j\in[0:N]}\lambda_{\star,j}\left\{ (\mathbf{g}_{\star}-\mathbf{g}_{j})\odot(\mathbf{g}_{\star}-\mathbf{g}_{j})\right\} \Bigg]\\
 & -\lambda_{\star,0}\{\mathbf{g}_{0}\odot\mathbf{w}_{0}\}\\
 & -\sum_{i\in[1:N-1]}\Bigg[\lambda_{i,i+1}\left\{ \mathbf{g}_{i}\odot\mathbf{w}_{0}\right\} -\sum_{j\in[0:i-1]}\frac{\alpha_{i,j}^{\prime}}{L}\left\{ \mathbf{g}_{i}\odot\mathbf{g}_{j}\right\} \Bigg]\\
 & -\left((\mathbf{g}_{N}\odot\mathbf{w}_{0})-\sum_{j\in[0:N-1]}\frac{\alpha_{N,j}^{\prime}}{L}\left\{ \mathbf{g}_{N}\odot\mathbf{g}_{j}\right\} \right)\\
 & +\sum_{i\in[0:N-1]}\Bigg[\lambda_{i,i+1}\left\{ \mathbf{g}_{i+1}\odot\mathbf{w}_{0}\right\} -\sum_{j\in[0:i-1]}\frac{\alpha_{i,j}^{\prime}}{L}\left\{ \mathbf{g}_{i+1}\odot\mathbf{g}_{j}\right\} \Bigg],
\end{align*}
$$
where

$$
\alpha_{i,j}^{\prime}=\begin{cases}
\lambda_{i,i+1}\alpha_{i,j}, & \textrm{if }i\in[1:N-1],j\in[0:i-1]\\
\alpha_{N,j}, & \textrm{if }i\in N,j\in[0:N-1]\\
\lambda_{0,1}\underbrace{\alpha_{0,j}}_{=0}=0, & \textrm{if }i=0.
\end{cases}
$$
Then, the performance estimation optimization algorithm is:
$$
\begin{array}{ll}
\textrm{minimize} & \tau\\
\textrm{subject to} & -\mathbf{f}_{N}+\mathbf{f}_{\star}+\sum_{(i,j)\in J}\lambda_{i,j}\left(\mathbf{f}_{j}-\mathbf{f}_{i}\right)+c_{f}\tau\left(\mathbf{f}_{0}-\mathbf{f}_{\star}\right)=0\\
 & S\left(\tau,\{\lambda_{i,j}\},\{\alpha_{i,j}^{\prime}\}\right)\succeq0\\
 & \left(\forall(i,j)\in J\right)\quad\lambda_{i,j}\geq0\\
 & \tau\geq0,
\end{array}
$$
where the decision variables are $\tau,\{\lambda_{i,j}\},$ and $\{\alpha_{i,j}^{\prime}\}.$



```julia
# Load the packages
using JuMP, MosekTools, Mosek, LinearAlgebra, SCS, COSMO, Literate, OffsetArrays
```

```julia
# %% Some helper functions

# %% construct e_i in R^n
function e_i(n, i)
    e_i_vec = zeros(n, 1)
    e_i_vec[i] = 1
    return e_i_vec
end

# %% construct symmetric outer product

function ⊙(a,b)
    return ((a*b')+(b*a')) ./ 2
end
```



```julia
# %% Parameters to be tuned
N = 5
L = 1
μ = 0
c_w = 1
c_f = 0
```

**Putting everything in a function.** Now, let us put everything in a function. We need to add the packages first.

```julia
function full_pep_solver(N, L, μ, c_w, c_f)

    # define all the bold vectors
    # --------------------------

    # define 𝐰_0

    𝐰_0 = e_i(N+2, 1)

    𝐰_star = zeros(N+2, 1)

    𝐠_star = zeros(N+2,1)

    # define 𝐠_0, 𝐠_1, …, 𝐠_N

    # first we define 𝐠_Julia vectors and then 𝐠 vectors

    𝐠 = OffsetArray(zeros(N+2, N+1), 1:N+2, 0:N)
    # 𝐠= [𝐠_0 𝐠_1 𝐠_2 ... 𝐠_N]
    for i in 0:N
        𝐠[:,i] = e_i(N+2, i+2)
    end


    𝐟_star = zeros(N+1,1)

    # time to define 𝐟_0, 𝐟_1, …, 𝐟_N

    𝐟 = OffsetArray(zeros(N+1, N+1), 1:N+1, 0:N)
    # 𝐟 = [𝐟_0, 𝐟_1, …, 𝐟_N]

    for i in 0:N
        𝐟[:,i] = e_i(N+1, i+1)
    end

    # Define JuMP model
    # ----------------

    pep_model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-10))

    #define all the variables
    # -------------------------------

    # define α′ (can be typed as \alpha[TAB]\prime[TAB])
    @variable(pep_model, α′[1:N, 0:N-1])

    # define λ variables

    @variable(pep_model, λ_i_ip1[0:N-1] >= 0)
    # this defines (λ_{i,i+1})_{i∈[0:N-1]} in Julia indexing

    @variable(pep_model, λ_star_i[0:N] >= 0)
    # this defines (λ_{⋆,i})_{i∈[0:N]} in Julia indexing

    # define τ
    @variable(pep_model, τ >= 0)

    # define objective
    # ------------------

    @objective(
    pep_model,
    Min,
    τ
    )

    # Add the linear equality constraint
    # ----------------------------------
    @constraint(pep_model,
    τ * c_f .* 𝐟[:,0] + sum(λ_i_ip1[i] .* (𝐟[:,i+1]-𝐟[:,i]) for i in 0:N-1)
    + sum(λ_star_i[i] .* (𝐟[:,i] - 𝐟_star) for i in 0:N)
    .== 𝐟[:,N] - 𝐟_star
    )

    # add the LMI constraint now

    @constraint(pep_model,
    # term_1:
    c_w * τ * ⊙(𝐰_0,𝐰_0) +
    # term_2 (part 1+part 2):
    (sum(λ_i_ip1[i] .* ⊙(𝐠[:,i]-𝐠[:,i+1],𝐠[:,i]-𝐠[:,i+1]) for i in 0:N-1) +
    sum(λ_star_i[i] .* ⊙(𝐠_star - 𝐠[:,i],𝐠_star - 𝐠[:,i]) for i in 0:N) ) ./ (2*L) +
    # term_3:
    (- λ_star_i[0] .* ⊙(𝐠[:,0],𝐰_0) ) +
    # term_4 part 1:
    (- sum( λ_i_ip1[i] .* ⊙(𝐠[:,i],𝐰_0) for i in 1:N-1)) +
    # term_4 part 2:
    (1/L) .* sum( sum(α′[i,j] .* ⊙(𝐠[:,i],𝐠[:,j]) for j in 0:i-1) for i in 1:N-1) +
    # term 5:
    (- (⊙(𝐠[:,N],𝐰_0) - (1/L) .* sum(α′[N,j] .* ⊙(𝐠[:,N],𝐠[:,j]) for j in 0:N-1))) +
    # term_6 part 1:
    sum(λ_i_ip1[i] .* ⊙(𝐠[:,i+1],𝐰_0) for i in 0:N-1) +
    # term_6 part 2:
    (- (1/L) .* sum( sum( α′[i,j] .* ⊙(𝐠[:,i+1],𝐠[:,j]) for j in 0:i-1) for i in 1:N-1))
    in PSDCone()
    )

    # time to optimize!
    # -----------------
    optimize!(pep_model)

    println("termination status =", termination_status(pep_model) )

    α′_opt = OffsetArray(value.(α′), 1:N, 0:N-1)
    λ_opt_i_ip1 = OffsetVector(value.(λ_i_ip1), 0:N-1)
    λ_opt_star_i = OffsetVector(value.(λ_star_i), 0:N)
    τ_opt = value(τ)

    # Recover α_{i,j} from α′_{i,j}
    # -----------------------------

    α_opt = OffsetArray(zeros(size(α′_opt)), 1:N, 0:N-1)

    for i in 1:N-1
        for j in 0:i-1
            α_opt[i,j] = (α′_opt[i,j] / λ_opt_i_ip1[i])
        end
    end

    for j in 0:N-1
        α_opt[N,j] = α′_opt[N,j]
    end

    # Recover h_{i,j} from α_{i,j}
    # ----------------------------

    h_opt =  OffsetArray(zeros(size(α_opt)), 1:N, 0:N-1)
    
    # set h(1,0)
    h_opt[1,0] = α_opt[1,0]

    # set the rest of the variables

    for i in 2:N
        for j in 0:i-1
            if j == i-1
                h_opt[i,j] = α_opt[i,j]
            else
                h_opt[i,j] = α_opt[i,j] - α_opt[i-1,j]
            end
        end
    end

    # return all the outputs
    # ----------------------

    return   α′_opt, λ_opt_i_ip1, λ_opt_star_i, τ_opt, α_opt, h_opt

end
```

Time to run and test.

```julia
α′_opt, λ_opt_i_ip1, λ_opt_star_i, τ_opt, α_opt, h_opt = full_pep_solver(N, L, μ, c_w, c_f)
```
