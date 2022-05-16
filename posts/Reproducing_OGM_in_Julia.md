@def title = "Reproducing Optimized Gradient Method (OGM) in Julia+JuMP"   
@def published = "May 16, 2022"   
@def tags =["programming", "Julia", "performance-estimation-problem"]  

# Reproducing Optimized Gradient Method (OGM) in Julia+JuMP

**Shuvomoy Das Gupta**

*May16, 2022*

In this blog, we discuss how to reproduce the Optimized Gradient Method (OGM) using `Julia+JuMP` step by step. OGM is the optimal method for minimizing smooth convex functions with respect to the performance measure of function value gap. It was first discovered numerically by Drori and Teboulle in their paper [https://arxiv.org/pdf/1206.3209.pdf](https://arxiv.org/pdf/1206.3209.pdf) and with the analytical form it appearing in the paper https://arxiv.org/pdf/1406.5468.pdf by Kim and Fessler. For more details, please see the aforementioned papers.

## Background

**Optimization problem.** The optimization problem is
$$
\begin{array}{ll}

\textrm{minimize} & f(x)\\

\textrm{subject to} & x\in\mathbf{R}^{d},

\end{array}
$$


where $f\in\mathcal{F}_{0,L}$ is the class of $L$-smooth convex functions.

**Fixed-step first-order algorithm.** The optimization algorithm in consideration is:
$$
\begin{align*}

\left(\forall_{i\in[0:N]}\right)\;w_{i} & =w_{i-1}-\sum_{j\in[0:i-1]}\frac{h_{i,j}}{L}\nabla f(w_{j})\\

 & =w_{0}-\sum_{i\in[0:i-1]}\frac{\alpha_{i,j}}{L}\nabla f(w_{j}),

\end{align*}
$$


where we use the notation that $\alpha_{0,i}=h_{0,i}=0$ and $\{\alpha_{i,j}\}$ and $\{h_{i,j}\}$ are related via the following relationship:


$$
\left(\forall i\in[1:N]\right)\;\left(\forall j\in[0:i-1]\right)\quad h_{i,j}=\begin{cases}

\alpha_{i,j}, & \textrm{if }j=i-1\\

\alpha_{i,j}-\alpha_{i-1,j}, & \textrm{if }j\in[0:i-1].

\end{cases}
$$



**Notation.** Denote,
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
where the decision variables are $\tau,\{\lambda_{i,j}\},$ and $\{\alpha_{i,j}^{\prime}\}.$  Let us see, how we can solve this problem in `Julia`.

## Reproducing OGM step-by-step in Julia+JuMP

```julia
# Load the packages
using JuMP, MosekTools, Mosek, LinearAlgebra, OffsetArrays
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
Next, define the bold vectors:

$$
\mathbf{w}_{0}=e_{1}\in\mathbf{R}^{N+2},\mathbf{g}_{i}=e_{i+2}\in\mathbf{R}^{N+2},\mathbf{f}_{i}=e_{i+1}\in\mathbf{R}^{N+1},
$$
where $e_{i}$ is the unit vector with $i$th component equal to $1$ and the rest being zero.

```julia
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
```
Next, we define our decision variables using `JuMP`.

```julia
pep_model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-10))

# time to define all the variables
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
```

Define the objective next, which is to minimize $\tau$.

```julia
# define objective
@objective(
    pep_model,
    Min,
    τ
    )
```

 Time to add the linear constraint
$$
\begin{align*}
 & -\mathbf{f}_{N}+\mathbf{f}_{\star}+\sum_{(i,j)\in J}\lambda_{i,j}\left(\mathbf{f}_{j}-\mathbf{f}_{i}\right)+c_{f}\tau\left(\mathbf{f}_{0}-\mathbf{f}_{\star}\right)=0\\
\Leftrightarrow & c_{f}\tau\left(\mathbf{f}_{0}-\mathbf{f}_{\star}\right)+\sum_{i\in[0:N-1]}\lambda_{i,i+1}(\mathbf{f}_{i+1}-\mathbf{f}_{i})+\sum_{i\in[0:N]}\lambda_{\star,i}(\mathbf{f}_{i}-\mathbf{f}_{\star})=\mathbf{f}_{N}-\mathbf{f}_{\star}
\end{align*}
$$
 first, which are much simpler.

```julia
# Add the linear equality constraint
# ----------------------------------
@constraint(pep_model,
    τ * c_f .* 𝐟[:,0] + sum(λ_i_ip1[i] .* (𝐟[:,i+1]-𝐟[:,i]) for i in 0:N-1)
    + sum(λ_star_i[i] .* (𝐟[:,i] - 𝐟_star) for i in 0:N)
    .== 𝐟[:,N] - 𝐟_star
)
```

Now let us construct the giant sdp constraint step by step. It has 6 summands, which we add one by one.

```julia
# Add the giant LMI constraint step by step
# ----------------------------------------
```
Term 1 is $\texttt{term}_1 = c_{w}\tau\mathbf{w}_{0}\mathbf{w}_{0}^{\top}.$
```julia
term_1 = c_w * τ * ⊙(𝐰_0,𝐰_0)
```
Term 2 is
$$
\texttt{term}_2 = \frac{1}{2L}\Bigg[\sum_{i\in[0:N-1]}\lambda_{i,i+1}\left\{ (\mathbf{g}_{i}-\mathbf{g}_{i+1})\odot(\mathbf{g}_{i}-\mathbf{g}_{i+1})\right\} +\\\sum_{j\in[0:N]}\lambda_{\star,j}\left\{ (\mathbf{g}_{\star}-\mathbf{g}_{j})\odot(\mathbf{g}_{\star}-\mathbf{g}_{j})\right\} \Bigg].
$$

```julia
term_2_part_1 = sum(λ_i_ip1[i] .* ⊙(𝐠[:,i]-𝐠[:,i+1],𝐠[:,i]-𝐠[:,i+1]) for i in 0:N-1)
```

```julia
term_2_part_2 = sum(λ_star_i[i] .* ⊙(𝐠_star - 𝐠[:,i],𝐠_star - 𝐠[:,i]) for i in 0:N)
```

```julia
term_2 = (term_2_part_1 + term_2_part_2) ./ (2*L)
```
Term 3 is $\texttt{term}_3 = -\lambda_{\star,0}\{\mathbf{g}_{0}\odot\mathbf{w}_{0}\}$.

```julia
term_3 = - λ_star_i[0] .* ⊙(𝐠[:,0],𝐰_0)
```
Term 4 is
$$
\texttt{term}_4 = -\sum_{i\in[1:N-1]}\Bigg[\lambda_{i,i+1}\left\{ \mathbf{g}_{i}\odot\mathbf{w}_{0}\right\} -\sum_{j\in[0:i-1]}\frac{\alpha_{i,j}^{\prime}}{L}\left\{ \mathbf{g}_{i}\odot\mathbf{g}_{j}\right\} \Bigg].
$$


```julia
term_4_part_1 = - sum( λ_i_ip1[i] .* ⊙(𝐠[:,i],𝐰_0) for i in 1:N-1)
```

```julia
term_4_part_2 = (1/L) .*
sum( sum(α′[i,j] .* ⊙(𝐠[:,i],𝐠[:,j]) for j in 0:i-1) for i in 1:N-1)
```

```julia
term_4 = term_4_part_1 + term_4_part_2
```
Term 5 is
$$
\texttt{term}_5 = -\left((\mathbf{g}_{N}\odot\mathbf{w}_{0})-\sum_{j\in[0:N-1]}\frac{\alpha_{N,j}^{\prime}}{L}\left\{ \mathbf{g}_{N}\odot\mathbf{g}_{j}\right\} \right).
$$


```julia
term_5 = - ( ⊙(𝐠[:,N],𝐰_0) - (1/L) .* sum(α′[N,j] .* ⊙(𝐠[:,N],𝐠[:,j]) for j in 0:N-1))
```
Okay, we are almost there. One more term. Term 6 is:
$$
\texttt{term}_6 = \sum_{i\in[0:N-1]}\Bigg[\lambda_{i,i+1}\left\{ \mathbf{g}_{i+1}\odot\mathbf{w}_{0}\right\} -\sum_{j\in[0:i-1]}\frac{\alpha_{i,j}^{\prime}}{L}\left\{ \mathbf{g}_{i+1}\odot\mathbf{g}_{j}\right\} \Bigg].
$$

```julia
term_6_part_1 = sum(λ_i_ip1[i] .* ⊙(𝐠[:,i+1],𝐰_0) for i in 0:N-1)
```

```julia
term_6_part_2 = - (1/L) .*
sum( sum( α′[i,j] .* ⊙(𝐠[:,i+1],𝐠[:,j]) for j in 0:i-1) for i in 1:N-1)
```

```julia
term_6 = term_6_part_1 + term_6_part_2
```

Time to add all these terms to construct $S\left(\tau,\{\lambda_{i,j}\},\{\alpha_{i,j}^{\prime}\}\right)$:
$$
S\left(\tau,\{\lambda_{i,j}\},\{\alpha_{i,j}^{\prime}\}\right) = \sum_{i\in [1:6]}\texttt{term}_i
$$

```julia
# oof, okay constructed the terms for the LMI constraint, let us hope that there is no mistake and add them together

S_mat = term_1 + term_2 + term_3 + term_4 + term_5 + term_6
```

Let's add the LMI cosntraint $S\left(\tau,\{\lambda_{i,j}\},\{\alpha_{i,j}^{\prime}\}\right)\succeq0$.

```julia
# add the LMI constraint now

@constraint(pep_model,
S_mat in PSDCone()
)
```

Time to optimize the code now!

```julia
# time to optimize!
# -----------------
optimize!(pep_model)

println("termination status =", termination_status(pep_model) )
```

The output is as follows.

```julia
# Problem
#   Name                   :
#   Objective sense        : min
#   Type                   : CONIC (conic optimization problem)
#   Constraints            : 34
#   Cones                  : 0
#   Scalar variables       : 37
#   Matrix variables       : 1
#   Integer variables      : 0
#
# Optimizer started.
# Presolve started.
# Linear dependency checker started.
# Linear dependency checker terminated.
# Eliminator started.
# Freed constraints in eliminator : 5
# Eliminator terminated.
# Eliminator started.
# Freed constraints in eliminator : 0
# Eliminator terminated.
# Eliminator - tries                  : 2                 time                   : 0.00
# Lin. dep.  - tries                  : 1                 time                   : 0.05
# Lin. dep.  - number                 : 0
# Presolve terminated. Time: 0.08
# Problem
#   Name                   :
#   Objective sense        : min
#   Type                   : CONIC (conic optimization problem)
#   Constraints            : 34
#   Cones                  : 0
#   Scalar variables       : 37
#   Matrix variables       : 1
#   Integer variables      : 0
#
# Optimizer  - threads                : 4
# Optimizer  - solved problem         : the primal
# Optimizer  - Constraints            : 29
# Optimizer  - Cones                  : 1
# Optimizer  - Scalar variables       : 23                conic                  : 16
# Optimizer  - Semi-definite variables: 1                 scalarized             : 28
# Factor     - setup time             : 0.02              dense det. time        : 0.00
# Factor     - ML order time          : 0.00              GP order time          : 0.00
# Factor     - nonzeros before factor : 423               after factor           : 435
# Factor     - dense dim.             : 0                 flops                  : 1.21e+04
# ITE PFEAS    DFEAS    GFEAS    PRSTATUS   POBJ              DOBJ              MU       TIME
# 0   1.0e+00  1.0e+00  1.0e+00  0.00e+00   0.000000000e+00   0.000000000e+00   1.0e+00  0.20
# 1   2.7e-01  2.7e-01  1.2e-01  2.14e+00   -2.378627944e-02  -6.117274954e-03  2.7e-01  0.28
# 2   1.2e-01  1.2e-01  2.9e-02  1.32e+00   7.175923608e-02   7.139302023e-02   1.2e-01  0.30
# 3   4.2e-02  4.2e-02  6.9e-03  1.05e+00   4.477758503e-02   4.729110878e-02   4.2e-02  0.30
# 4   2.3e-02  2.3e-02  2.9e-03  7.37e-01   2.880293993e-02   3.057913998e-02   2.3e-02  0.31
# 5   6.7e-03  6.7e-03  4.0e-04  9.85e-01   3.356646000e-02   3.302278094e-02   6.7e-03  0.33
# 6   2.5e-03  2.5e-03  1.1e-04  6.92e-01   2.291866424e-02   2.288414204e-02   2.5e-03  0.33
# 7   6.6e-04  6.6e-04  1.5e-05  8.71e-01   2.066303721e-02   2.066465445e-02   6.6e-04  0.33
# 8   9.8e-05  9.8e-05  8.8e-07  8.83e-01   1.883068878e-02   1.882786667e-02   9.8e-05  0.33
# 9   1.7e-06  1.7e-06  1.9e-09  9.93e-01   1.859217137e-02   1.859212284e-02   1.7e-06  0.33
# 10  4.0e-08  4.0e-08  7.3e-12  1.00e+00   1.858820828e-02   1.858820710e-02   4.0e-08  0.34
# 11  2.1e-09  2.1e-09  8.5e-14  1.00e+00   1.858813924e-02   1.858813917e-02   2.1e-09  0.34
# 12  7.9e-11  1.8e-10  6.3e-16  1.00e+00   1.858813673e-02   1.858813673e-02   7.9e-11  0.34
# Optimizer terminated. Time: 0.45
#
# termination status =OPTIMAL
```

We now collect the optimal values of our decision variables.

```julia
# Collect the decision variables
# -----------------------------
α′_opt = OffsetArray(value.(α′), 1:N, 0:N-1)
λ_opt_i_ip1 = OffsetVector(value.(λ_i_ip1), 0:N-1)
λ_opt_star_i = OffsetVector(value.(λ_star_i), 0:N)
τ_opt = value(τ)
```

The resultant  is as follows.

```julia
# α′_opt =
# [
# 0.314962  0.0       0.0      0.0      0.0
# 0.641151  0.722442  0.0      0.0      0.0
# 1.05006   1.38407   1.2547   0.0      0.0
# 1.54002   2.17685   2.32945  1.90951  0.0
# 1.92565   2.8008    3.17533  2.96989  2.07777
# ]
```



Next, we recover $\alpha_{i,j}$ from $\alpha^\prime_{i,j}$ using the formula
$$
\alpha_{i,j}=\begin{cases}
\frac{\alpha_{i,j}^{\prime}}{\lambda_{i,i+1}}, & \textrm{if }i\in[1:N-1],j\in[0:i-1]\\
\alpha_{N,j}^{\prime}, & \textrm{if }i=N.
\end{cases}
$$

```julia
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
```
The output for `α_opt` is as follows.

```julia
# α_opt
# =
# [
#  1.61803  0.0      0.0      0.0      0.0
#  1.79217  2.01939  0.0      0.0      0.0
#  1.86775  2.46185  2.23175  0.0      0.0
#  1.90789  2.69683  2.88589  2.36563  0.0
#  1.92565  2.8008   3.17533  2.96989  2.07777
#  ]
```

We now, recover $h_{i,j}$ from $\alpha_{i,j}$ using the formula
$$
\left(\forall i\in[1:N]\right)\;\left(\forall j\in[0:i-1]\right)\quad h_{i,j}=\begin{cases}
\alpha_{i,j}, & \textrm{if }j=i-1\\
\alpha_{i,j}-\alpha_{i-1,j}, & \textrm{if }j\in[0:i-1].
\end{cases}
$$

```julia
# Recover h_{i,j} from α_{i,j}
h_opt =  OffsetArray(zeros(size(α_opt)), 1:N, 0:N-1)

# set h(1,0)
h_opt[1,0] = α_opt[1,0]

for i in 2:N
    for j in 0:i-1
        if j == i-1
            h_opt[i,j] = α_opt[i,j]
        else
            h_opt[i,j] = α_opt[i,j] - α_opt[i-1,j]
        end
    end
end
```

We have the following values for `h_opt`.

```julia
# h_opt
# =
# [
# 1.61803    0.0       0.0       0.0       0.0
# 0.174133   2.01939   0.0       0.0       0.0
# 0.0755813  0.442461  2.23175   0.0       0.0
# 0.0401385  0.234975  0.654138  2.36563   0.0
# 0.0177604  0.103971  0.289442  0.604262  2.07777
# ]
```

## Putting everything in a compact function

Now, let us put everything in a function. We need to add the packages first.

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

    # Add the giant LMI constraint step by step
    # ----------------------------------------

    # Define all the terms one by one

    term_1 = c_w * τ * ⊙(𝐰_0,𝐰_0)

    term_2_part_1 = sum(λ_i_ip1[i] .* ⊙(𝐠[:,i]-𝐠[:,i+1],𝐠[:,i]-𝐠[:,i+1]) for i in 0:N-1)

    term_2_part_2 = sum(λ_star_i[i] .* ⊙(𝐠_star - 𝐠[:,i],𝐠_star - 𝐠[:,i]) for i in 0:N)

    term_2 = (term_2_part_1 + term_2_part_2) ./ (2*L)

    term_3 = - λ_star_i[0] .* ⊙(𝐠[:,0],𝐰_0)

    term_4_part_1 = - sum( λ_i_ip1[i] .* ⊙(𝐠[:,i],𝐰_0) for i in 1:N-1)

    term_4_part_2 = (1/L) .*
    sum( sum(α′[i,j] .* ⊙(𝐠[:,i],𝐠[:,j]) for j in 0:i-1) for i in 1:N-1)

    term_4 = term_4_part_1 + term_4_part_2

    term_5 = - ( ⊙(𝐠[:,N],𝐰_0) - (1/L) .* sum(α′[N,j] .* ⊙(𝐠[:,N],𝐠[:,j]) for j in 0:N-1))

    term_6_part_1 = sum(λ_i_ip1[i] .* ⊙(𝐠[:,i+1],𝐰_0) for i in 0:N-1)

    term_6_part_2 = - (1/L) .*
    sum( sum( α′[i,j] .* ⊙(𝐠[:,i+1],𝐠[:,j]) for j in 0:i-1) for i in 1:N-1)

    term_6 = term_6_part_1 + term_6_part_2

    # oof, okay constructed the terms for the LMI constraint, let us hope that there is no mistake and add them together

    S_mat = term_1 + term_2 + term_3 + term_4 + term_5 + term_6

    # add the LMI constraint now

    @constraint(pep_model,
    S_mat in PSDCone()
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
