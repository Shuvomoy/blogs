@def title = "Computing proximal operator of a constrained function in Julia"
@def published ="September 24, 2020"
@def tags =["programming", "Julia"]

# Computing proximal operator of a constrained function in Julia

**Shuvomoy Das Gupta**

*September 24, 2020*

---

**Table of contents**
\toc

---

(Last update: January 13, 2021) In this blog, we will show how to compute proximal operator of a constrained function. The entire code is written in markdown using the package [Weave.jl](https://github.com/JunoLab/Weave.jl). The corresponding `.jmd` notebook (can be opened both as a markdown file in any text editor or code in `Juno`), can be downloaded [here](https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/codes/2020-09-08-proximal_operator_over_matrix.jmd).

We will consider two examples, the first one arises in low-rank factor analysis problem, and the second one comes from lifted problem of a standard form binary integer optimization problem.

### Example 1: Low-rank factor analysis

As an example we consider the function:

$$
f(X,D) =\left\Vert \Sigma-X-D\right\Vert _{F}^{2}+I_{\mathcal{P}}(X,D),
$$

where $I_{\mathcal{P}}$ denotes the indicator function of the convex
set
$$
\mathcal{P}=\{(X,D)\in\mathbf{S}^{n}\times\mathbf{S}^{n}\mid X\succeq0,D=\mathbf{diag}(d),d\geq0,d\in \mathbf{R}^{n}, \Sigma-D \succeq 0\}.
$$

#### Computing the proximal operator of $f$

Proximal operator $\mathbf{prox}_{\gamma f}$ for this function $f$ at $(X,D)$ is *the* optimal solution to the following convex optimization problem:

$$
\begin{equation}
\begin{array}{ll}
\textrm{minimize} & \left\Vert \Sigma-\widetilde{X}-\widetilde{D}\right\Vert _{F}^{2}+\frac{1}{2\gamma}\|\widetilde{X}-X\|_{F}^{2}+\frac{1}{2\gamma}\|\widetilde{D}-D\|_{F}^{2}\\
\textrm{subject to} & \widetilde{X}\succeq0\\
 & \widetilde{D}=\mathbf{diag}(\widetilde{d})\\
 & \widetilde{d}\geq0, \\
 & \Sigma - \widetilde{D} \succeq 0
\end{array}
\end{equation}
$$

where $\widetilde{X}\in\mathbf{S}_{+}^{n},$ and $\widetilde{d}\in \mathbf{R}_{+}^{n}$(*i.e.*, $\widetilde{D}=\mathbf{diag}(\widetilde{d}$)) are the optimization variables.

Now we solve this optimization problem using `Julia`. We will use the package `Convex` and `SCS`, both open source `Julia` packages.

#### Load the packages

First, we load the packages. If the packages are not installed we can install them by running the following commands in `Julia` REPL.

```julia
using Pkg
Pkg.add("Convex")
Pkg.add("COSMO")
```


```julia
## Load the packages
using Convex
using LinearAlgebra
using COSMO
using JuMP
using SCS
```

#### Solver function

Let us write the solver function that is going to solve the optimization problem that we described above. The first implementation is using `Convex.jl` and the second one is via `JuMP`.

##### First implementation using `Convex.jl`


```julia
# put everything in a function: first implementatino is using Convex.jl
function prox_PRS_fam_cvxjl(Σ, M, γ, X, d) #(Σ::A, M::R, γ::R, X::A, d::V) where {R <: Real, A <: AbstractMatrix{R}, V <:  AbstractVector{R}} # For now M is not used, may use it in a future version

  # This functions takes the input data Σ, γ, X, d and evaluates
  # the proximal operator of the function f at the point (X,D)

  # Data extraction
  # ---------------
  n = length(d) # dimension of the problem
  # println("*****************************")
  # println(size(d))
  # println("the value of d is = ", d)
  # println("the type of d is", typeof(d))
  D = LinearAlgebra.diagm(d) # creates the diagonal matrix D that embed

  # Create the variables
  #  --------------------
  X_tl = Convex.Semidefinite(n) # Here Semidefinite(n) encodes that
  # X_tl ≡ ̃X is a positive semidefinite matrix
  d_tl = Convex.Variable(n) # d_tl ≡ ̃d
  D_tl = diagm(d_tl) # Create the diagonal matrix ̃D from ̃d

  # Create terms of the objective function, which we write down
  #  in three parts
  #  ----------------------------------------------------------
  t1 = square(norm2(Σ - X_tl - D_tl)) # norm2 computes Frobenius norm in Convex.jl
  t2 = square(norm2(X-X_tl))
  t3 = square(norm2(D-D_tl))

  # Create objective
  # ----------------
  objective = t1 + (1/(2*γ))*(t2 + t3) # the objective to be minimized

  # create the problem instance
  # ---------------------------
  convex_problem = Convex.minimize(objective, [d_tl >= 0, Σ - D_tl  in :SDP])

  # set the solver
  # --------------
  convex_solver = () -> SCS.Optimizer(verbose=false)

  # solve the problem
  # -----------------
  Convex.solve!(convex_problem, convex_solver)

  # get the optimal solution
  # ------------------------
  X_sol = X_tl.value
  d_sol = d_tl.value
  # println("d_sol = ", d_sol)

  # return the output
  return X_sol, vec(d_sol)

end
```

##### Second implementation using `JuMP`

```julia
## put everything in a function (implementation using JuMP)
function prox_PRS_fam_JuMP(Σ, M, γ, X, d; X_tl_sv = nothing, d_tl_sv = nothing, warm_start = false)

	# This functions takes the input data Σ, γ, X, d and evaluates the proximal operator of the function f at the point (X,d)

	n = length(d)

	prox_model = JuMP.Model(optimizer_with_attributes(SCS.Optimizer, "verbose" => false))

	# prox_model = JuMP.Model(optimizer_with_attributes(Mosek.Optimizer))


	@variables( prox_model,
	begin
		d_tl[1:n] >= 0
		X_tl[1:n, 1:n], PSD
	end
	)

	if warm_start == true
		println("warm start enabled")
		# Set warm-start
		set_start_value.(X_tl, X_tl_sv) # Warm start
		set_start_value.(d_tl, d_tl_sv) # Warm start
	    println("norm difference is = ", norm(start_value.(X_tl) - X_tl_sv))
	end



    t_1 = vec(Σ - X_tl - diagm(d_tl))
	t_2 = vec(X_tl-X)
	t_3 = vec(diagm(d_tl)-diagm(d))
	obj = t_1'*t_1 + ((1/(2*γ))*(t_2'*t_2 + t_3'*t_3))

	@objective(prox_model, Min, obj)

	@constraints(prox_model, begin
		Symmetric(Σ - diagm(d_tl)) in PSDCone()
	end)

	set_silent(prox_model)

	JuMP.optimize!(prox_model)

	# obj_val = JuMP.objective_value(prox_model)
	X_sol = JuMP.value.(X_tl)
	d_sol = JuMP.value.(d_tl)

	return X_sol, d_sol

end
```




#### Create data to test

We create the data now to test the function we just wrote.


```julia
## All the parameters
n = 10
Σ1 = randn(n,n)
Σ = Σ1'*Σ1
X = randn(n,n)
d = randn(n)
M = 1
γ = 1
```

#### Test the function

We test the function now to see if the function `prox_over_matrix` works as expected!


```julia
## Time to run the code
@time X1, d1 = prox_PRS_fam_cvxjl(Σ, M, γ, X, d)

@time X2, d2 = prox_PRS_fam_JuMP(Σ, M, γ, X, d)

## Test for warmstarting

X_tl_sv = X2
d_tl_sv = d2

@time X2, d2 = prox_PRS_fam_JuMP(Σ, M, γ, X, d; X_tl_sv = X2, d_tl_sv = d2, warm_start = true)
```

`Output:`

```julia
  2.263578 seconds (11.69 k allocations: 1.033 MiB)
  1.613454 seconds (18.60 k allocations: 2.731 MiB)
  warm start enabled
  norm difference is = 0.0
  1.643022 seconds (19.30 k allocations: 2.764 MiB)
```

### Example 2. Lifted SDP relaxation of a standard form 0-1 integer optimization problem

The standard form binary integer optimization problem
is a universal combinatorial problem with the following form:
$$
\begin{equation}
\begin{array}{ll}
\textrm{minimize} & c^{\top}x\\
\textrm{subject to} & Ax=b\\
 & x\geq0\\
 & x\in\{0,1\}^{n},
\end{array}
\end{equation}
$$
where $x$ is the decision variable, and $A\in\mathbf{R}^{m\times n},b\in\mathbf{R}^{m},c\in\mathbf{R}^{n}$ are problem data. The matrix $A$ is taken to be fat $(m<n)$ and full row rank. This problem has the following higher-dimensional equivalent formulation based on the lifting procedure originally proposed by
Lovász and Schrijver:
$$
\begin{equation}
\begin{array}{ll}
\textrm{minimize} & \sum_{j=1}^{n}c_{j}X_{jj}\\
\textrm{subject to} & \sum_{j=1}^{n}A_{j}X_{ij}=bX_{ii},\qquad i=1,\ldots,n,\\
 & \sum_{j=1}^{n}A_{j}X_{jj}=b\\
 & 0\leq X_{ij}\leq X_{ii},\qquad i,j=1,\ldots,n,i\neq j,\\
 & 0\leq X_{ij}\leq1,\qquad i,j=1,\ldots,n,\\
 & X_{ii}+X_{jj}-X_{ij}\leq1,\qquad i,j=1,\ldots,n,i\neq j,\\
 & X\succeq0,\\
 & \mathbf{rank} X=1,
\end{array}
\end{equation}
$$

where $X\in\mathbf{S}^{n}$ is the decision variable and is related to the original decision variable $x$ by the encoding $X=xx^{\top}.$ In this formulation $A_{j}$ represents the $j$-th column of the matrix $A$.

Let us drop the nonconvex rank constraint, and consider the function
$$
f(X) :=\sum_{j=1}^{n}c_{j}X_{jj} +I_{\mathcal{C}}(X),
$$

where $I_{\mathcal{C}}$ denotes the indicator function of the convex set

$$
\begin{equation}
\begin{array}{ll}
\mathcal{C} := & \{ X\in\mathbf{S}^{n}\mid\sum_{j=1}^{n}A_{j}X_{ij}=bX_{ii},\qquad i\in1,\ldots,n,\\
 & \qquad\qquad\sum_{j=1}^{n}A_{j}X_{jj}=b,\\
 & \qquad\qquad 0\leq X_{ij}\leq X_{ii},\qquad i,j=1,\ldots,n,i\neq j,\\
 & \qquad \qquad0\leq X_{ij}\leq1,\qquad i,j=1,\ldots,n\\
 & \qquad \qquad X_{ii}+X_{jj}-X_{ij}\leq1,\qquad i,j=1,\ldots,n,i\neq j,\\
 & \qquad \qquad X\succeq0 \}.
\end{array}
\end{equation}
$$

#### Computing the proximal operator of $f$

Proximal operator $\mathbf{prox}_{\gamma f}$ for this function $f$ at $X$ is *the* optimal solution to the following convex optimization problem:

$$
\begin{equation}
\begin{array}{ll}
\textrm{minimize} & \sum_{j=1}^{n}c_{j}\widetilde{X}_{jj}+\frac{1}{2\gamma}\|\widetilde{X}-X\|_{F}^{2}\\
\textrm{subject to} & \sum_{j=1}^{n}A_{j}\widetilde X_{ij}=b \widetilde X_{ii},\qquad i=1,\ldots,n,\\
 & \sum_{j=1}^{n}A_{j} \widetilde X_{jj}=b\\
 & 0\leq \widetilde X_{ij}\leq \widetilde X_{ii},\qquad i,j=1,\ldots,n,i\neq j,\\
 & 0\leq \widetilde X_{ij}\leq1,\qquad i,j=1,\ldots,n,\\
 & \widetilde X_{ii}+\widetilde X_{jj}-\widetilde X_{ij}\leq1,\qquad i,j=1,\ldots,n,i\neq j,\\
 & \widetilde X\succeq0,\\

\end{array}
\end{equation}
$$

where $\widetilde{X}\in\mathbf{S}_{+}^{n}$ is the optimization variables.

#### Load the packages

First, we load the packages.

```julia
## Load the packages
using Convex
using LinearAlgebra
using COSMO
using JuMP
using SCS
```

```julia
# Problem Data
m = 5
n = 10
using Random: bitrand
c = randn(n)
x0 = bitrand(n)
A = randn(m,n)
b = A*x0 # this is to ensure that the problem has a solution
x_test = randn(n)
X = x_test*x_test' # this is the point where we want to evaluate the proximal operator
γ = 1 # proximal parameter
```

```julia
# Evaluating the proximal operator (put the code in a function)
prox_model = JuMP.Model(optimizer_with_attributes(Mosek.Optimizer))

@variables(prox_model,
begin
	X_tl[1:n, 1:n], PSD
end
)

t_1 = sum(c[j]*X_tl[j,j] for j in 1:n)
t_2 = vec(X_tl-X)
obj = t_1 + ((1/(2*γ))*(t_2'*t_2))

@objective(prox_model, Min, obj)

@constraints(prox_model, begin
con_sdp_equality_1[i = 1:n], sum(X_tl[j,j].*A[:,j] for j in 1:n) .== b*X_tl[i,i]
con_sdp_equality_2, sum(X_tl[j,j].*A[:,j] for j in 1:n) .== b
con_bound[i = 1:n, j = 1:n], 0<= X_tl[i,j] <= 1
cons_sdp_bound_1[i = 1:n, j= 1:n, i != j], X_tl[i,i]+X_tl[j,j]-X_tl[i,j] <= 1
con_sdp_bound_2[i = 1:n, j = 1:n, i != j], X_tl[i,j] <= X_tl[i,i]
end)

set_silent(prox_model)

JuMP.optimize!(prox_model)

obj_val = JuMP.objective_value(prox_model)

X_sol = JuMP.value.(X_tl)
```

###### Converting `.jmd` file to `.jl` file

To convert the `.jmd` file to a `.jl` file we run the following code:

```julia
using Weave
cd("C:\\Users\\shuvo\\Google Drive\\GitHub\\blog\\codes") # directory that contains the .jmd file
tangle("2020-09-08-proximal_operator_over_matrix.jmd", informat = "markdown") # convert the .jmd file into a .jl file that will contain the code
```