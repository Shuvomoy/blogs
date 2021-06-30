@def title = "  Solving mixed-integer quadratic programming (MIQP) using Julia+JuMP"
@def published ="January 20, 2020"
@def tags =["programming", "Julia", "optimization"]

In this blog, we will discuss how to solve a mixed-integer quadratic programming problem (MIQP) using `Julia` and `JuMP`. My versions of `Julia`, `JuMP`, and `Gurobi` are 1.3.0, 0.20.1, and 0.7.4, respectively.

As an illustrative example, we will consider the sparse regression problem. The sparse regression is a nonconvex optimization problem with applications to gene expression analysis and signal processing etc.

**Sparse regression problem.** The sparse regression problem is concerned with approximating a vector $b\in\mathbf{R}^{m}$ with a linear combination of at most $k$ columns of a matrix $A\in\mathbf{R}^{m\times n}$ with bounded coefficients. The problem can be written as the following optimization problem  
$$
\begin{array}{ll} \textrm{minimize} & \|Ax-b\|^{2}\\ \textrm{subject to} & \mathbf{card}(x)\leq k\\ & \|x\|_{\infty}\leq M, \end{array} 
$$


where $x\in\mathbf{R}^{n}$ is the decision variable, and $A\in\mathbf{R}^{m\times n},b\in\mathbf{R}^{m},$ and $M\in\mathbf{R}$ are problem data. Here, $\mathbf{card}(x)$ is the number of nonzero components in $x$.

**Modeling sparse regression problem as a MIQP.** The sparse regression problem can be modeled as the following MIQP: 
$$
\begin{array}{ll} \textrm{minimize} & \|Ax-b\|^{2}\\ \textrm{subject to} & |x_{i}|\leq My_{i},\quad i=1,\ldots,n \qquad (1)\\ & \sum_{i=1}^{n}y_{i}\leq k\\ & x\in\mathbf{R}^{n},y\in\{0,1\}^{n}, \end{array}
$$
 where $x, y$ are decision variables.

We can write our objective function as 


$$
\begin{align*} \|Ax-b\|^{2} & =(Ax-b)^{\intercal}(Ax-b)\\ & =(x^{\intercal}A^{\intercal}-b^{\intercal})(Ax-b)\\ & =x^{\intercal}(A^{\intercal}A)x+(-2A^{\intercal}b)^{\intercal}x+\|b\|^{2}, 
\end{align*}
$$
which takes the function in a quadratic form.

**Rewriting the objective in a compatible format.** If we define $S=A^{\intercal}A,$ $c=-2A^{\intercal}b,$ and $d=\|b\|^{2},$ then we can write the objective function as: $$ \sum_{i=1}^{n}\sum_{j=1}^{n}S_{ij}x_{i}x_{j}+\sum_{i=1}^{n}c_{i}x_{i}+d, $$ which is more compatible as an input for `JuMP`.

**Rewriting the bound constraint in a more compatible format.** Also, for `JuMP`, we write the bound constraint $|x_{i}|\leq My_{i}$ as two constraints: $x_i \leq M y_i$, and $-M y_i \leq x_i$ for $i=1,\ldots,n$.

The code is as follows.

In [1]:

```
# Load the necessary packages
# ---------------------------

using Gurobi, JuMP, LinearAlgebra
# if these packages are not installed then run:
# using Pkg
# Pkg.add("Gurobi")
# Pkg.add("JuMP")
# Also, keep in mind that Gurobi is a commercial solver, but it is free for academic use.
```

In [2]:

```
# Data, change it accordingly
# ---------------------------

m = 5
n = 10
A = randn(m,n)
b = randn(m)
M = 1
k = convert(Int64, round(m/3))

# Renaming a bunch of variables
S = A'*A
c = -2*A'*b
d = norm(b)^2
```

Out[2]:

```
3.511217774252138
```

In [3]:

```
# Define the model
# ----------------

model = Model(with_optimizer(Gurobi.Optimizer)) # define name of the model, it could be anything, not necessarily "model"

# Variables
# ---------

@variable(model, x[1:n]) # define variable x 

@variable(model, y[1:n], Bin) # define the binary variable y

# Objective
# ---------

sense = MOI.MIN_SENSE # by this command, we are programatically defining a quadratic objective to be minimized 

@objective(model, sense, sum(S[i,j]*x[i]*x[j] for i in 1:n, j in 1:n)+ sum(c[i]*x[i] for i in 1:n) + d) # define the objective

# Constraints
# ------------

@constraint(model, con_lb[i=1:n], -M*y[i] <= x[i]) # lower bound constraint

@constraint(model, con_ub[i=1:n], x[i] <= M*y[i]) # upper bound constraint

@constraint(model, con_bd_sum, sum(y[i] for i in 1:n) <= k) # cardinality constraint in terms of y

# Run the optimizer
# -----------------

status=optimize!(model) # time to optimize!

# Let us look at the important outputs 
# ------------------------------------
println("******************************************************")
println("optimal objective value is = ", objective_value(model))
println("optimal x is = ",  value.(x))
println("optimal y is =", value.(y))
Academic license - for non-commercial use only
Academic license - for non-commercial use only
Optimize a model with 21 rows, 20 columns and 50 nonzeros
Model has 55 quadratic objective terms
Variable types: 10 continuous, 10 integer (10 binary)
Coefficient statistics:
  Matrix range     [1e+00, 1e+00]
  Objective range  [6e-01, 7e+00]
  QObjective range [2e-02, 3e+01]
  Bounds range     [0e+00, 0e+00]
  RHS range        [2e+00, 2e+00]
Found heuristic solution: objective 3.5112178
Presolve time: 0.00s
Presolved: 21 rows, 20 columns, 50 nonzeros
Presolved model has 55 quadratic objective terms
Variable types: 10 continuous, 10 integer (10 binary)

Root relaxation: objective 1.776357e-15, 45 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0    0.00000    0    6    3.51122    0.00000   100%     -    0s
H    0     0                       1.5146229    0.00000   100%     -    0s
     0     0    0.00000    0    6    1.51462    0.00000   100%     -    0s
H    0     0                       0.9576732    0.00000   100%     -    0s
     0     2    0.00000    0    6    0.95767    0.00000   100%     -    0s
H   18     7                       0.7568978    0.59570  21.3%   4.8    0s

Explored 25 nodes (161 simplex iterations) in 0.02 seconds
Thread count was 8 (of 8 available processors)

Solution count 4: 0.756898 0.957673 1.51462 3.51122 

Optimal solution found (tolerance 1.00e-04)
Best objective 7.568977587066e-01, best bound 7.568977587066e-01, gap 0.0000%
******************************************************
optimal objective value is = 0.7568977587066126
optimal x is = [-1.0, 0.0, -0.9733324670061584, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
optimal y is =[1.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0]
```