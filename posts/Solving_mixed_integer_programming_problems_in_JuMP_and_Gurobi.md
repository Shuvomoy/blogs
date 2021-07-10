@def title = "Solving mixed integer programming (MIP) problems using Julia+JuMP"
@def published ="July 10, 2021"
@def tags =["programming", "Julia", "optimization"]



# Solving mixed integer programming (MIP) problems using Julia+JuMP

**Shuvomoy Das Gupta**

*July 10, 2021*

In this blog, we will discuss how to solve a mixed integer programming (MIP) problem using Julia and JuMP.

---

\toc

---

### Problem

Let us try to write the JuMP code for the following standard form optimization problem:

$$
\begin{align}
& \text{minimize} && c^T x + d^T y\\
& \text{subject to} && A x + B y= f \\
 &                   && x \succeq 0, y \succeq 0 \\
 &                   && x \in \mathbb{R}^n, y \in \mathbb{\{0,1\}}^p,
\end{align}
$$

where $x,y$ are the decision variables, $A \in \mathbb{R}^{m \times n}, B \in \mathbb{R}^{m \times p}, c \in \mathbb{R}^n, d \in \mathbb{R}^p, f \in \mathbb{R}^m$. The data were randomly generated. The symbol $\succeq$ ($\preceq$) stands for element-wise greater (less) than or equal to.

### Data generation

Let us generate random data.

```julia
using Random

n = 15
p = 14
m = 13

A = randn(m,n)

B = randn(m,p)

c = abs.(randn(n,1))

d = abs.(randn(p,1))

x_rand_feas = abs.(randn(n,1))

y_rand_feas = bitrand(p,1)

f = A*x_rand + B*y_rand # to ensure that we have a feasible solution

```
```julia
using JuMP, Gurobi
```

### Function for solving MIP

Let us write the function, which will find an optimal solution. We have used the notions of warm-start and variable hint in the function. Information about them can be found at the following links:

* [Warm-start](https://www.gurobi.com/documentation/9.1/refman/start.html),

* [Variable hint](https://www.gurobi.com/documentation/9.1/refman/varhintval.html#attr:VarHintVal).

```julia
function MIP_solver(A, B, c, d, f;
    # if warm-start value is provided
    warm_start = :off,
    x_start_value = zeros(n, 1),
    y_start_value = zeros(p,1),
    # if variable hint is provided
    variable_hint = :off,
    x_hint =  zeros(n, 1),
    y_hint = zeros(p,1)
    )

    sfMipModel = direct_model(Gurobi.Optimizer())
    # using direct_model has results in smaller memory allocation


    set_optimizer_attribute(sfMipModel, "MIPFocus", 3)
    # If you are more interested in good quality feasible solutions, you can select MIPFocus=1.
    # If you believe the solver is having no trouble finding the optimal solution, and wish to focus more
    # attention on proving optimality, select MIPFocus=2.
    # If the best objective bound is moving very slowly (or not at all), you may want to try MIPFocus=3 to focus on the bound.


    set_optimizer_attribute(sfMipModel, "NonConvex", 2)
    # "NonConvex" => 2 tells Gurobi to use its nonconvex algorithm

    # some other Gurobi options relevant for MIP solving (commented here)
    # -------------------------------------------------------------------

    # set_optimizer_attribute(BnB_PEP_model, "MIPGap", 1e-4)
    ## for more info see
    ## https://www.gurobi.com/documentation/9.1/refman/mipgap2.html#parameter:MIPGap

    # set_optimizer_attribute(BnB_PEP_model, "FeasibilityTol", 1e-4)
    ## for more info see
    ## https://www.gurobi.com/documentation/9.1/refman/feasibilitytol.html

    # set_optimizer_attribute(BnB_PEP_model, "OptimalityTol", 1e-4)
    ## for more info see
    ## https://www.gurobi.com/documentation/9.1/refman/optimalitytol.html

    @variable(sfMipModel, x[1:n] >=0)
    @variable(sfMipModel, y[1:p] >= 0, Bin)

    if warm_start == :on

        for i in 1:n
            set_start_value(x[i], x_start_value[i])
            # alternatively the following code will do the same thing
            # MOI.set(sfMipModel, Gurobi.VariableAttribute("Start"), x[i], x_start_value[i])
        end

        for i in 1:p
            set_start_value(y[i], y_start_value[i])
            # alternatively the following code will do the same thing
            # MOI.set(sfMipModel, Gurobi.VariableAttribute("Start"), y[i], y_start_value[i])
        end

    end

    if variable_hint == :on
        for i in 1:n
            MOI.set(sfMipModel, Gurobi.VariableAttribute("VarHintVal"), x[i], x_start_value[i])
            MOI.set(sfMipModel, Gurobi.VariableAttribute("VarHintPri"), x[i], 10) # this number represents our confidence in the variable hint
        end

        for i in 1:p
            MOI.set(sfMipModel, Gurobi.VariableAttribute("VarHintVal"), y[i], y_start_value[i])
            MOI.set(sfMipModel, Gurobi.VariableAttribute("VarHintPri"), y[i], 10) # this number represents our confidence in the variable hint
        end

    end

    @objective(sfMipModel, Min, sum(c[i] * x[i] for i in 1:n)+sum(d[i]*y[i] for i in 1:p))

    for i in 1:m
        @constraint(sfMipModel, sum(A[i,j]*x[j] for j in 1:n)+ sum(B[i,j]*y[j] for j in 1:p) == f[i])
    end

    optimize!(sfMipModel)

    # store the important values and return them
    status = JuMP.termination_status(sfMipModel)
    x_sol = value.(x)
    y_sol = value.(y)
    obj_value = objective_value(sfMipModel)
    
    return obj_value, x_sol, y_sol

end
```

### Solving from scratch

```julia
obj_value_1, x_sol_1, y_sol_1 = MIP_solver(A, B, c, d, f)
```

```julia
Solver output
-------------
Gurobi Optimizer version 9.1.2 build v9.1.2rc0 (win64)
Thread count: 4 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 13 rows, 29 columns and 377 nonzeros
Model fingerprint: 0xf6c02cdc
Variable types: 15 continuous, 14 integer (14 binary)
Coefficient statistics:
  Matrix range     [2e-03, 3e+00]
  Objective range  [4e-03, 2e+00]
  Bounds range     [0e+00, 0e+00]
  RHS range        [6e-01, 1e+01]
Presolve time: 0.00s
Presolved: 13 rows, 29 columns, 356 nonzeros
Variable types: 15 continuous, 14 integer (14 binary)
Presolved: 13 rows, 29 columns, 356 nonzeros


Root relaxation: objective 4.410336e+00, 16 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0    4.41034    0    5          -    4.41034      -     -    0s
H    0     0                       7.0242168    4.41034  37.2%     -    0s
H    0     0                       6.5680788    4.41034  32.9%     -    0s
     0     0    6.56808    0    4    6.56808    6.56808  0.00%     -    0s

Cutting planes:
  Gomory: 5
  MIR: 3
  Flow cover: 6

Explored 1 nodes (24 simplex iterations) in 0.06 seconds
Thread count was 8 (of 8 available processors)

Solution count 2: 6.56808 7.02422

Optimal solution found (tolerance 1.00e-04)
Best objective 6.568078761353e+00, best bound 6.568078761353e+00, gap 0.0000%

User-callback calls 74, time in user-callback 0.00 sec
```

### Solving using warm-start

```julia
obj_value_2, x_sol_2, y_sol_2 = MIP_solver(A, B, c, d, f;
    # if warm-start value is provided
    warm_start = :on,
    x_start_value = x_sol_1,
    y_start_value = y_sol_1 )
```
```julia
Gurobi Optimizer version 9.1.2 build v9.1.2rc0 (win64)
Thread count: 4 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 13 rows, 29 columns and 377 nonzeros
Model fingerprint: 0xbeda5cf1
Variable types: 15 continuous, 14 integer (14 binary)
Coefficient statistics:
  Matrix range     [2e-03, 3e+00]
  Objective range  [4e-03, 2e+00]
  Bounds range     [0e+00, 0e+00]
  RHS range        [6e-01, 1e+01]

Loaded user MIP start with objective 6.56808 # this tells us that our provided warm-start values has been accepted by Gurobi as a feasible solution

Presolve time: 0.00s
Presolved: 13 rows, 29 columns, 356 nonzeros
Variable types: 15 continuous, 14 integer (14 binary)
Presolved: 13 rows, 29 columns, 356 nonzeros


Root relaxation: objective 4.410336e+00, 16 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0    4.41034    0    5    6.56808    4.41034  32.9%     -    0s

Cutting planes:
  Gomory: 5
  MIR: 3
  Flow cover: 6

Explored 1 nodes (16 simplex iterations) in 0.08 seconds
Thread count was 8 (of 8 available processors)

Solution count 1: 6.56808

Optimal solution found (tolerance 1.00e-04)
Best objective 6.568078761353e+00, best bound 6.568078761353e+00, gap 0.0000%

User-callback calls 71, time in user-callback 0.00 sec
```

### Solving using variable hint

```julia
obj_value_3, x_sol_3, y_sol_3 = MIP_solver(A, B, c, d, f;
    # if variable hint is provided
    variable_hint = :on,
    x_hint =  x_sol_1,
    y_hint = y_sol_1)
```

```julia
Solver output
-------------
Gurobi Optimizer version 9.1.2 build v9.1.2rc0 (win64)
Thread count: 4 physical cores, 8 logical processors, using up to 8 threads
Optimize a model with 13 rows, 29 columns and 377 nonzeros
Model fingerprint: 0x9ff923d5
Variable types: 15 continuous, 14 integer (14 binary)
Coefficient statistics:
  Matrix range     [2e-03, 3e+00]
  Objective range  [4e-03, 2e+00]
  Bounds range     [0e+00, 0e+00]
  RHS range        [6e-01, 1e+01]
Presolve time: 0.00s
Presolved: 13 rows, 29 columns, 356 nonzeros
Variable types: 15 continuous, 14 integer (14 binary)
Presolved: 13 rows, 29 columns, 356 nonzeros


Root relaxation: objective 4.410336e+00, 16 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0    4.41034    0    5          -    4.41034      -     -    0s
H    0     0                       7.0242168    4.41034  37.2%     -    0s
H    0     0                       6.5680788    4.41034  32.9%     -    0s
     0     0    6.56808    0    4    6.56808    6.56808  0.00%     -    0s

Cutting planes:
  Gomory: 5
  MIR: 3
  Flow cover: 6

Explored 1 nodes (24 simplex iterations) in 0.05 seconds
Thread count was 8 (of 8 available processors)

Solution count 2: 6.56808 7.02422

Optimal solution found (tolerance 1.00e-04)
Best objective 6.568078761353e+00, best bound 6.568078761353e+00, gap 0.0000%

User-callback calls 71, time in user-callback 0.00 sec
```


### Some helpful online links on speeding up the solution process

1. [How do I use MIP starts](https://support.gurobi.com/hc/en-us/articles/360043834831-How-do-I-use-MIP-starts-)
2. [Warmstart or VarHintVal in QP Python API](https://support.gurobi.com/hc/en-us/community/posts/360075269952-Warmstart-or-VarHintVal-in-QP-Python-API)
3. [Most important parameters for MIP models](https://www.gurobi.com/documentation/9.1/refman/mip_models.html)
