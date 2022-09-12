@def title = "Solving bilinear optimization problem using Alpine.jl"
@def published ="August 30, 2022"
@def tags =["programming", "Julia", "optimization"]

# Solving bilinear optimization problem using Alpine.jl

**Shuvomoy Das Gupta**

*August 30, 2022*

In this blog, we will discuss how to solve bilinear optimizations problem using `Alpine`. 

---

**Table of contents**

\toc

---

## What is `Alpine`?

`Alpine` is an open-source `Julia` package to solve **m**ixed-**i**nteger **n**on**l**inear optimization **p**roblems (MINLPs) to global optimality. It is an open-source implementation of the [**A**daptive **M**ultivariate **P**artitioning (AMP) Algorithm](https://arxiv.org/pdf/1707.02514.pdf) proposed by Nagarajan et al in [Nagarajan2019]. AMP uses an adaptive, piecewise convexification scheme and constraint programming methods to solve MINLPs to global optimality. The benefit of `Alpine` over other spatial branch-and-bound solvers is that it is entirely built upon `JuMP` [Dunning2017] and `MathOptInterface` [Legat2021] in Julia, which provides significant flexibility and modeling power. 

## Example in consideration

We want to solve the following bilinear problem:
$$
\begin{array}{ll}
\textrm{maximize} & x_{1}x_d - x_2 x_{d-1}\\
\textrm{subject to} & \sum_{i=1}^{d}x_{i}\leq10\\
 & x_{i}x_{i+1}\leq2,\quad i=1,\ldots,d-1\\
 & \sum_{i=1}^{d-1}x_{i}x_{i+1}=1\\
 & x\geq0,
\end{array} \quad \textrm{(BLP)}
$$
where $x\in\mathbf{R}^d$ is the decision variables. For this example, we let $d=10$

Now, we will solve (BLP) step by step in `Alpine`.

## Construct the optimization model using `JuMP`

First, we load the packages.

```julia
## Load packages
using Alpine, JuMP, Gurobi, Ipopt
```

Now, we create a function that will construct the optimization model in  `JuMP`. For the optimization model to be classified as an MINLP, we need at least one nonlinear constraint or objective declared specifically in the `JuMP` model. This can be achieved by using `@NLobjective` or `@NLconstraint` macro.

```julia
function BLPSolver(; solver = AlpineSolver) 
# solver is either AlpineSolver or GurobiSolver (defined below)

    d = 10
  
    # declare the JuMP model
    # ----------------------
    m = Model(solver)

    # declare the variable x ≧ 0
    # ----------------------
    @variable(m, x[1:d] >= 0) 
  
    # declare the constraints
    # -----------------------

    # add objective x[1] x[d] - x[2] x[d-1], which is to be maximized

    if solver ==  GurobiSolver
        @objective(m, Max, x[1] * x[d] - x[2] * x[d-1])
    elseif solver == AlpineSolver
        @NLobjective(m, Max, x[1] * x[d] - x[2] * x[d-1])
    else 
       @error "correct solver type not given"
    end

    # add the constraints
    # -------------------

    # constraint ∑_{i=1:d} x[i] ≦ 10, conSum is the name of the constraint
    @constraint(m, conSum, sum(x[i] for i in 1:d) <= 10) 

    # add the constraints  (∀i ∈ [1:d-1]) x[i] x[i+1] ≤ 2, 
    # where the i-the constraint is called congBilinear1[i]
    @constraint(m, conBilinear1[i=1:d-1], x[i] * x[i+1] <= 2)

    # add constraint ∑_{i ∈ [1:d-1]} x[i] x[i+1] = 1, we call it 
    @constraint(m, conBilinear2, sum(x[i] * x[i+1] for i in 1:d-1) == 1)
  
    # return the JuMP model
    # ---------------------

    return m

end
```

## Setting local and global subsolvers in `Alpine`

Now we will solve this problem using `Alpine`. To implement the MAP algorithm, which is the core algorithm in  `Alpine` , we need:

(i) One local solver, which usually implements some variant of the nonlinear interior point algorithm to find a locally optimal solution to the MINLP. Solvers of this type that are supported by `Alpine` via ` MathOptInterface` are: `Ipopt` and `KNITRO`.

(ii) One mixed-integer optimization solver, that usually implements some sort of branch-and-bound algorithm to solve mixed-integer optimization problems (MIPs). During every iteration of the MAP algorithm, `Alpine` solves an MIP sub-problem. Because MIPs are usually hard to solve, `Alpine`'s runtime heavily depends on the runtime of these solvers. The MIP solvers that are supported by `Alpine` are: `Cbc`, `CPLEX`, `Gurobi`,  and `Bonmin`. In my numerical experiments, I have found `Gurobi` to be the fastest and the `Alpine`  developers also recommend `Gurobi` as well.

We set the local and global subsolvers, followed by setting `Alpine` as the global solver as follows.  

```julia
# function to set local subsolver as Ipopt
# ----------------------------------------
IpoptSolverAlpine = JuMP.optimizer_with_attributes(
                                      Ipopt.Optimizer,
                                      MOI.Silent() => true
                                      )

# function to set global subsolver as Gurobi
# ------------------------------------------
GRB_ENV  = Gurobi.Env()

GurobiSolverAlpine = JuMP.optimizer_with_attributes(
                                        () -> Gurobi.Optimizer(GRB_ENV),
                                        MOI.Silent() => true,
                                        "Presolve" => 0
                                        )

# function to set Alpine as the global solver
# -------------------------------------------
AlpineSolver = JuMP.optimizer_with_attributes(
                                              Alpine.Optimizer,
                                              "nlp_solver"   => IpoptSolverAlpine, # local solver
                                              "mip_solver"   => GurobiSolverAlpine, # global solver
                                              "presolve_bt"  => true, 
                                              # "disc_ratio"   => 10,
                                              "apply_partitioning" => true,
                                              "log_level" => 1
                                              )                                    
```

Define a Gurobi solver to solve the bilinear problem to global optimality. 

```julia
GurobiSolver = JuMP.optimizer_with_attributes(
                                        () -> Gurobi.Optimizer(GRB_ENV),
                                        # MOI.Silent() => true,
                                        "NonConvex" => 2 
                                        # means we are going to use Gurobi's spatial branch-and-bound algorithm
                                        )
```

The options used in declaring the `AlpineSolver` solver above are as follows (Source: [https://lanl-ansi.github.io/Alpine.jl/latest/parameters/](https://lanl-ansi.github.io/Alpine.jl/latest/parameters/))

`"nlp_solver"   => IpoptSolverAlpine`  sets `Ipopt` as the local solver 

`"mip_solver"   => GurobiSolverAlpine`  sets `Gurobi` as the global solver

`"presolve_bt"  => true`  performs sequential optimization-based bound tightening (OBBT) at the presolve step. For more details about sequential OBBT, see Section 3.1.1 of [https://arxiv.org/pdf/1707.02514.pdf](https://arxiv.org/pdf/1707.02514.pdf).

`disc_ratio => 10`  is used to measure the width of new partitions relative to the active partition chosen in the sequentially solved lower-bounding MIPs. This value can substantially affect the run time for global convergence; this value can be set to different integer values (>= 4) for various classes of problems.

`"apply_partitioning" => true`  applies `Alpine`'s built-in MIP-based partitioning algorithm  only (MAP) when activated; else terminates with the presolve solution.

`"log_level" => 100`  enables detailed debugging mode of `Alpine`. The option `log_level (default = 0)` controls the verbosity level of Alpine output; choose 1 for turning on logging, else 100 for detailed debugging mode.

## Solve the problem using `Alpine`

Now we solve the problem using Alpine.

```julia
m = BLPSolver(solver = AlpineSolver)

JuMP.optimize!(m)
```

The output is as follows.

```
PROBLEM STATISTICS
  Objective sense = Max
  # Variables = 10
  # Bin-Int Variables = 0
  # Constraints = 11
  # NL Constraints = 10
  # Linear Constraints = 1
  # Detected convex constraints = 0
  # Detected nonlinear terms = 11
  # Variables involved in nonlinear terms = 10
  # Potential variables for partitioning = 10
SUB-SOLVERS USED BY ALPINE
  NLP local solver = Ipopt
  MIP solver = Gurobi
ALPINE CONFIGURATION
  Maximum iterations (upper-bounding MIPs) =  99
  Relative global optimality gap = 0.01%
  Discretization ratio = 10
  Bound-tightening presolve = true
  Maximum iterations (OBBT) = 25
PRESOLVE 
  Doing local search
  Local solver returns a feasible point with value 23.9792
  Starting bound-tightening
  Actual iterations (OBBT): 9
  Post-presolve optimality gap: 0.569%
  Completed presolve in 0.33s
UPPER-BOUNDING ITERATIONS
====================================================================================================
| Iter   | Incumbent       | Best Incumbent      | Upper Bound        | Gap (%)         | Time      
| 1      | 23.9792         | 23.9792             | 24.009             | 0.124           | 0.34s            
| 2      | 24.0            | 24.0                | 24.009             | 0.038           | 0.38s            
| finish | 24.0            | 24.0                | 24.0021            | 0.009           | 0.46s            
====================================================================================================

## *** Alpine ended with status OPTIMAL *** ##
```

## Solving the problem using `Gurobi`

Alternatively, we can also solve the problem in consideration to global optimality using `Gurobi` as follows. 


Now solve the problem. 

```julia
mGurobi = BLPSolver(solver = GurobiSolver);

JuMP.optimize!(mGurobi)
```

Output is as follows.

```
Presolve time: 0.00s
Presolved: 61 rows, 22 columns, 120 nonzeros
Presolved model has 20 bilinear constraint(s)
Variable types: 22 continuous, 0 integer (0 binary)

Root relaxation: objective 4.955556e+01, 22 iterations, 0.00 seconds (0.00 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   49.55556    0   11          -   49.55556      -     -    0s
H    0     0                      23.9791576   49.55556   107%     -    0s
     0     0   24.99322    0   10   23.97916   24.99322  4.23%     -    0s
     0     0   24.99231    0    8   23.97916   24.99231  4.23%     -    0s
     0     0   24.99029    0   10   23.97916   24.99029  4.22%     -    0s
     0     0   24.99026    0   11   23.97916   24.99026  4.22%     -    0s
     0     2   24.99026    0   11   23.97916   24.99026  4.22%     -    0s
H    5     6                      23.9946925   24.95311  3.99%  11.0    0s
*   37    10               7      23.9961128   24.05999  0.27%   6.5    0s
*   45    12               8      23.9972011   24.01419  0.07%   5.6    0s
*   50    12               9      23.9980299   24.01286  0.06%   5.3    0s
*   57    12               9      23.9998195   24.01286  0.05%   4.8    0s
H   64    14                      23.9999797   24.01166  0.05%   4.3    0s
H  141    65                      23.9999994   24.00776  0.03%   2.2    0s
*  172    47              24      24.0000004   24.00776  0.03%   1.9    0s
*  185    47              24      24.0000004   24.00776  0.03%   1.7    0s

Cutting planes:
  RLT: 13

Explored 562 nodes (382 simplex iterations) in 0.05 seconds (0.01 work units)
Thread count was 10 (of 10 available processors)

Solution count 10: 24 24 24 ... 23.9792

Optimal solution found (tolerance 1.00e-04)
Best objective 2.399999999991e+01, best bound 2.400042876397e+01, gap 0.0018%

User-callback calls 1347, time in user-callback 0.00 sec
```

Note that we have the same optimal value. 

## Getting tightened bounds on the variables using Alpine

In certain applications, we may be interested in computing valid but tightened bounds on the variables rather than finding an optimal solution. One nice feature of Alpine is that, just by using its sequential OBBT algorithm, we can compute valid tightened bounds on the variables without solving the full problem. We can achieve this in two ways. 



* The first way is to set `"apply_partitioning" => false`, which will terminate with the presolve solution but with tighter bounds. 



* The second way is keeping `"apply_partitioning" => true,` but setting the option `"rel_gap" => δ`, where δ is a number larger than or equal to the the default relative gap `1e-4`. This will lead to Alpine terminating at the relative gap δ, so with a suboptimal solution with optimality gap δ%, but with significantly tighter bounds on the variables. 

For example, by inspecting BLP, we see that $0 \leq x \leq 10$. We will now show how to tighten these bounds significantly by using the sequential OBBT feature of the AMP algorithm that is implemented in `Alpine`. 

### The first way

```julia
## Set Alpine as the global solver but to compute tighter bounds via its presolve phase
AlpineSolver = JuMP.optimizer_with_attributes(
                                               Alpine.Optimizer,
                                              "nlp_solver"   => IpoptSolverAlpine, # local solver
                                              "mip_solver"   => GurobiSolverAlpine, # global solver
                                              "presolve_bt"  => true, 
                                              "partition_scaling_factor"  => 10,
                                              "apply_partitioning" => false,
                                              "log_level" => 1
                                              )           
```

Solve the bound tightening problem and get the tightened bounds on the variables.

```julia
mOBBT1 = BLPSolver(solver = AlpineSolver)

JuMP.optimize!(mOBBT1)
```

Let us get the lower bound and upper bound on our decision variable `x`.

```julia
vars = all_variables(mOBBT1) # 

lowerBoundVar = MOI.get.(mOBBT1, Alpine.TightenedLowerBound(), vars) 

upperBoundVar = MOI.get.(mOBBT1, Alpine.TightenedUpperBound(), vars)

@info "[☼ ] Bound on the decision variable x is:"

for i in eachindex(vars)
    @info "$(lowerBoundVar[i]) <=  x[$(i)] <= $(upperBoundVar[i])"
end
```

The output is:

```
[ Info: [☼ ] Bound on the decision variable x is:

[ Info: 4.6326 <=  x[1] <= 5.174300000000001
[ Info: 0.0 <=  x[2] <= 0.2096
[ Info: 0.0 <=  x[3] <= 0.030600000000000002
[ Info: 0.0 <=  x[4] <= 0.030000000000000002
[ Info: 0.0 <=  x[5] <= 0.030000000000000002
[ Info: 0.0 <=  x[6] <= 0.030000000000000002
[ Info: 0.0 <=  x[7] <= 0.030000000000000002
[ Info: 0.0 <=  x[8] <= 0.030600000000000002
[ Info: 0.0 <=  x[9] <= 0.2096
[ Info: 4.633100000000001 <=  x[10] <= 5.1742
```

, which is significantly tighter than $0 \leq x \leq 10$.

### The second way 

```julia
## Set Alpine as the global solver but to compute tighter bounds
AlpineSolver = JuMP.optimizer_with_attributes(
                                               Alpine.Optimizer,
                                              "nlp_solver"   => IpoptSolverAlpine, # local solver
                                              "mip_solver"   => GurobiSolverAlpine, # global solver
                                              "presolve_bt"  => true, 
                                              "partition_scaling_factor"   => 10,
                                              "apply_partitioning" => true,
                                              "rel_gap" => 5e-2, # it means that we will stop when we have relative gap [(UpperBound-LowerBound)/UpperBound] <= 50%, because our objective is just getting tightened bounds on the variables rather than finding globally optimal solution
                                              "log_level" => 1
                                              )           
```

Solve the bound tightening problem and get the tightened bounds on the variables.

```julia
mOBBT2 = BLPSolver(solver = AlpineSolver)

JuMP.optimize!(mOBBT2)
```

Let us get the lower bound and upper bound on our decision variable `x`.

```julia
vars = all_variables(mOBBT2) # 

lowerBoundVar = MOI.get.(mOBBT2, Alpine.TightenedLowerBound(), vars) 

upperBoundVar = MOI.get.(mOBBT2, Alpine.TightenedUpperBound(), vars)

@info "[☼ ] Bound on the decision variable x is:"

for i in eachindex(vars)
    @info "$(lowerBoundVar[i]) <=  x[$(i)] <= $(upperBoundVar[i])"
end
```

The output is:

```
[ Info: [☼ ] Bound on the decision variable x is:
[ Info: 3.8099000000000003 <=  x[1] <= 5.7202
[ Info: 0.0 <=  x[2] <= 0.2625
[ Info: 0.0 <=  x[3] <= 0.334
[ Info: 0.0 <=  x[4] <= 0.329
[ Info: 0.0 <=  x[5] <= 0.329
[ Info: 0.0 <=  x[6] <= 0.329
[ Info: 0.0 <=  x[7] <= 0.329
[ Info: 0.0 <=  x[8] <= 0.3508
[ Info: 0.0 <=  x[9] <= 0.264
[ Info: 4.341 <=  x[10] <= 5.7131
```

, which is significantly tighter than $0 \leq x \leq 10$.



## References

[Nagarajan2019] Nagarajan, Harsha, et al. "An adaptive, multivariate partitioning algorithm for global optimization of nonconvex programs." Journal of Global Optimization 74.4 (2019): 639-675.

[Legat2021] Legat, Benoît, et al. "MathOptInterface: a data structure for mathematical optimization problems." INFORMS Journal on Computing 34.2 (2022): 672-689.

[Dunning2017] Dunning, Iain, Joey Huchette, and Miles Lubin. "JuMP: A modeling language for mathematical optimization." SIAM review 59.2 (2017): 295-320.

