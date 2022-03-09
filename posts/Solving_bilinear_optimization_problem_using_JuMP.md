@def title = "Solving bilinear optimization problem using JuMP"
@def published ="March 09, 2022"
@def tags =["programming", "Julia", "optimization"]

# Solving bilinear optimization problem using JuMP

**Shuvomoy Das Gupta**

*March 9, 2022*

In this blog, we will discuss how to solve a bilinear optimization problem  in `JuMP+Julia` using `Gurobi` and /or`KNITRO` as the solver. Here we note that, `Gurobi` is able to solve these problems to global optimality, whereas `KNITRO` can provide only solve these problems only locally. 

---

**Table of contents**

\toc

---

### Problem in consideration

We want to solve the following bilinear problem:
$$
\begin{array}{ll}
\textrm{maximize} & x_{1}x_d - x_2 x_{d-1}\\
\textrm{subject to} & \sum_{i=1}^{d}x_{i}\leq10\\
 & x_{i}x_{i+1}\leq2,\quad i=1,\ldots,d-1\\
 & \sum_{i=1}^{d-1}x_{i}x_{i+1}=1\\
 & x\geq0,
\end{array}
$$
where $x\in\mathbf{R}^d$ is the decision variables.

### `JuMP` code

The `JuMP`  code is as follows.

```julia
using Gurobi, JuMP, BenchmarkTools

function solve_bilinear(d; start_value_provided = :off, x_start_value = zeros(d), warm_start_type = :none, silent_flag = :off)

  # declare the model
  # -----------------

  nonlinear_model = direct_model(Gurobi.Optimizer())
  # using direct_model has results in smaller memory allocation
  # we could also use
  # nonlinear_model = Model(Gurobi.Optimizer)
  # but this requires more memory allocation

  set_optimizer_attribute(nonlinear_model, "MIPFocus", 3)
  # If you are more interested in good quality feasible solutions, you can select MIPFocus=1.
  # If you believe the solver is having no trouble finding the optimal solution, and wish to focus more
  # attention on proving optimality, select MIPFocus=2.
  # If the best objective bound is moving very slowly (or not at all), you may want to try MIPFocus=3 to focus on the bound.


  set_optimizer_attribute(nonlinear_model, "NonConvex", 2)
  # "NonConvex" => 2 tells Gurobi to use its nonconvex algorithm

  # declare the variable along with the non-negativity constraitn
  @variable(nonlinear_model, x[1:d] >=0 )

  # warm starting using Gurobi.jl
  if  start_value_provided == :on && warm_start_type == :Gurobi_direct
    for i in 1:d
        MOI.set(nonlinear_model, Gurobi.VariableAttribute("Start"), x[i], x_start_value[i])
    end
  end

  # warm starting using JuMP
  if  start_value_provided == :on && warm_start_type == :JuMP
    for i in 1:d
        set_start_value(x[i], x_start_value[i])
    end
  end

  # add objective
  @objective(nonlinear_model, Max, x[1]*x[d] - x[2]*x[d-1])

  # add the constraint Σx[i] <= 10
  @constraint(nonlinear_model, sum(x[i] for i in 1:d) <= 10)

  # add the bilinear inequality constraint ∀i  x[i]x[i+1] <= 2
  for i in 1:d-1
    @constraint(nonlinear_model, x[i]*x[i+1] <= 2)
  end

  # add the bilinear equality constraint ∑_{i} x[i]x[i+1] == 1
  @constraint(nonlinear_model, sum(x[i]*x[i+1] for i in 1:d-1) == 1)

  # optimize the model
  if silent_flag == :on
     set_silent(nonlinear_model)
   end
  optimize!(nonlinear_model)

  # store the important values and return them
  status = JuMP.termination_status(nonlinear_model)
  x_sol = JuMP.value.(x)
  obj_value = JuMP.objective_value(nonlinear_model)

  return status, x_sol, obj_value

end

```

### Solve from scratch

```julia
d = 10
status, x_sol, obj_value =  solve_bilinear(d; start_value_provided = :off, x_start_value = zeros(d), warm_start_type = :none, silent_flag = :off)
```
The output is:

```julia
Continuous model is non-convex -- solving as a MIP.

Presolve time: 0.02s
Presolved: 61 rows, 22 columns, 120 nonzeros
Presolved model has 20 bilinear constraint(s)
Variable types: 22 continuous, 0 integer (0 binary)
Presolve removed 39 rows and 2 columns
Presolved: 22 rows, 20 columns, 59 nonzeros


Root relaxation: objective 4.955556e+01, 19 iterations, 0.01 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   49.55556    0   11          -   49.55556      -     -    0s
     0     0   25.21339    0    2          -   25.21339      -     -    0s
     0     0   25.10018    0   11          -   25.10018      -     -    0s
     0     0   25.09269    0   11          -   25.09269      -     -    0s
     0     0   25.03708    0   10          -   25.03708      -     -    0s
     0     0   25.00988    0   10          -   25.00988      -     -    0s
     0     0   24.99742    0   10          -   24.99742      -     -    0s
     0     0   24.99291    0   10          -   24.99291      -     -    0s
H    0     0                    -100.0000000   24.99291   125%     -    0s
H    0     0                      23.9993796   24.99291  4.14%     -    0s
     0     2   24.99291    0   10   23.99938   24.99291  4.14%     -    0s

Cutting planes:
  RLT: 22

Explored 81 nodes (380 simplex iterations) in 0.18 seconds
Thread count was 8 (of 8 available processors)

Solution count 2: 23.9994 -100

Optimal solution found (tolerance 1.00e-04)
Best objective 2.399937963333e+01, best bound 2.400151700222e+01, gap 0.0089%

User-callback calls 298, time in user-callback 0.01 sec
```

Let us benchmark the results.

```julia
b1 = @benchmark solve_bilinear(d; start_value_provided = :off, x_start_value = zeros(d), warm_start_type = :none,  silent_flag = :on)
println("benchmark for original code without start value")
println("***************************")
io = IOBuffer()
show(io, "text/plain", b1)
s = String(take!(io))
println(s)
```

Benchmarked output is as follows.

```julia
benchmark for original code without start value
***************************
BenchmarkTools.Trial:
  memory estimate:  48.24 KiB
  allocs estimate:  884
  --------------
  minimum time:     45.491 ms (0.00% GC)
  median time:      55.348 ms (0.00% GC)
  mean time:        56.950 ms (0.00% GC)
  maximum time:     70.355 ms (0.00% GC)
  --------------
  samples:          88
  evals/sample:     1
```

### Warm-starting `Gurobi` with a known feasible solution

What if we have a good guess of a solution for this problem? In that case, we can feed that vaue using `set_start_value(variable_name, variable_start_value)` function or `MOI.set(nonlinear_model, Gurobi.VariableAttribute("Start"), variable_name, variable_start_value)`. To test this, we just feed the optimal solution that we just found to test if the warm-starting is working.

### Warm-starting via `JuMP`

```julia
status_2, x_sol_2, obj_value_2 = solve_bilinear(d; start_value_provided = :on, x_start_value = x_sol, warm_start_type = :JuMP, silent_flag = :off)
```

The output is:

```julia
Continuous model is non-convex -- solving as a MIP.


Loaded user MIP start with objective 23.9994

Presolve time: 0.00s
Presolved: 61 rows, 22 columns, 120 nonzeros
Presolved model has 20 bilinear constraint(s)
Variable types: 22 continuous, 0 integer (0 binary)
Presolve removed 39 rows and 2 columns
Presolved: 22 rows, 20 columns, 59 nonzeros


Root relaxation: objective 4.955556e+01, 19 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   49.55556    0   11   23.99938   49.55556   106%     -    0s
     0     0   25.21341    0    2   23.99938   25.21341  5.06%     -    0s
     0     0   25.09295    0    8   23.99938   25.09295  4.56%     -    0s
     0     0   25.03730    0   10   23.99938   25.03730  4.32%     -    0s
     0     0   25.01008    0   11   23.99938   25.01008  4.21%     -    0s
     0     0   24.99757    0    6   23.99938   24.99757  4.16%     -    0s
     0     0   24.99757    0   10   23.99938   24.99757  4.16%     -    0s
     0     0   24.99299    0    7   23.99938   24.99299  4.14%     -    0s
     0     2   24.99299    0    7   23.99938   24.99299  4.14%     -    0s

Cutting planes:
  RLT: 22

Explored 81 nodes (440 simplex iterations) in 0.08 seconds
Thread count was 8 (of 8 available processors)

Solution count 1: 23.9994

Optimal solution found (tolerance 1.00e-04)
Best objective 2.399937963333e+01, best bound 2.400120366930e+01, gap 0.0076%

User-callback calls 274, time in user-callback 0.00 sec
```

```julia
b2 = @benchmark solve_bilinear(d; start_value_provided = :on, x_start_value = x_sol, warm_start_type = :JuMP, silent_flag = :on)
println("benchmark for original code with start value via JuMP")
println("***************************")
io = IOBuffer()
show(io, "text/plain", b2)
s = String(take!(io))
println(s)
```

The benchmarking result is as follows.

```julia
benchmark for original code with start value via JuMP
***************************
BenchmarkTools.Trial:
  memory estimate:  48.43 KiB
  allocs estimate:  904
  --------------
  minimum time:     43.148 ms (0.00% GC)
  median time:      47.913 ms (0.00% GC)
  mean time:        49.426 ms (0.00% GC)
  maximum time:     84.259 ms (0.00% GC)
  --------------
  samples:          102
  evals/sample:     1
```

### Warm-starting via `Gurobi.jl`

```julia
status_3, x_sol_3, obj_value_3 = solve_bilinear(d; start_value_provided = :on, x_start_value = x_sol, warm_start_type =  :Gurobi_direct, silent_flag = :off)
```

The output is:

```julia
Continuous model is non-convex -- solving as a MIP.


Loaded user MIP start with objective 23.9994

Presolve time: 0.00s
Presolved: 61 rows, 22 columns, 120 nonzeros
Presolved model has 20 bilinear constraint(s)
Variable types: 22 continuous, 0 integer (0 binary)
Presolve removed 39 rows and 2 columns
Presolved: 22 rows, 20 columns, 59 nonzeros


Root relaxation: objective 4.955556e+01, 19 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   49.55556    0   11   23.99938   49.55556   106%     -    0s
     0     0   25.21341    0    2   23.99938   25.21341  5.06%     -    0s
     0     0   25.09295    0    8   23.99938   25.09295  4.56%     -    0s
     0     0   25.03730    0   10   23.99938   25.03730  4.32%     -    0s
     0     0   25.01008    0   11   23.99938   25.01008  4.21%     -    0s
     0     0   24.99757    0    6   23.99938   24.99757  4.16%     -    0s
     0     0   24.99757    0   10   23.99938   24.99757  4.16%     -    0s
     0     0   24.99299    0    7   23.99938   24.99299  4.14%     -    0s
     0     2   24.99299    0    7   23.99938   24.99299  4.14%     -    0s

Cutting planes:
  RLT: 22

Explored 81 nodes (440 simplex iterations) in 0.09 seconds
Thread count was 8 (of 8 available processors)

Solution count 1: 23.9994

Optimal solution found (tolerance 1.00e-04)
Best objective 2.399937963333e+01, best bound 2.400120366930e+01, gap 0.0076%

User-callback calls 267, time in user-callback 0.00 sec
```

Let's take a look at the benchmarking output for this case.

```julia
b3 = @benchmark solve_bilinear(d; start_value_provided = :on, x_start_value = x_sol, warm_start_type = :Gurobi_direct, silent_flag = :on)
println("benchmark for original code with start value via Gurobi.jl")
println("***************************")
io = IOBuffer()
show(io, "text/plain", b3)
s = String(take!(io))
println(s)
```

The output is as follows this time.

```julia
benchmark for original code with start value via Gurobi.jl
***************************
BenchmarkTools.Trial:
  memory estimate:  48.46 KiB
  allocs estimate:  910
  --------------
  minimum time:     31.733 ms (0.00% GC)
  median time:      36.254 ms (0.00% GC)
  mean time:        37.185 ms (0.00% GC)
  maximum time:     73.067 ms (0.00% GC)
  --------------
  samples:          135
  evals/sample:     1
```
We see that in both warm-starting procedures, `Gurobi` accepts the provided solution, which makes the branch and bound method work somewhat better (number of explored nodes goes down from 397 to 225 after warm-starting). Of course, because the problem size is small, it does not impact the runtime significantly.

### Finding a locally optimal solution using `KNITRO`

```julia
using JuMP, KNITRO
function solve_bilinear_locally_KNITRO(d; start_value_provided = :off, x_start_value = zeros(d), warm_start_type = :none, silent_flag = :off)

  # declare the model
  # -----------------

  nonlinear_model = Model(optimizer_with_attributes(KNITRO.Optimizer, "convex" => 0))

  # declare the variable along with the non-negativity constraitn
  @variable(nonlinear_model, x[1:d] >=0 )


  # warm starting using JuMP
  if  start_value_provided == :on
    for i in 1:d
        set_start_value(x[i], x_start_value[i])
    end
  end

  # add objective
  @objective(nonlinear_model, Max, x[1]*x[d] - x[2]*x[d-1])

  # add the constraint Σx[i] <= 10
  @constraint(nonlinear_model, sum(x[i] for i in 1:d) <= 10)

  # add the bilinear inequality constraint ∀i  x[i]x[i+1] <= 2
  for i in 1:d-1
    @constraint(nonlinear_model, x[i]*x[i+1] <= 2)
  end

  # add the bilinear equality constraint ∑_{i} x[i]x[i+1] == 1
  @constraint(nonlinear_model, sum(x[i]*x[i+1] for i in 1:d-1) == 1)

  # optimize the model
  if silent_flag == :on
     set_silent(nonlinear_model)
   end
  optimize!(nonlinear_model)

  # store the important values and return them
  status = JuMP.termination_status(nonlinear_model)
  x_sol = JuMP.value.(x)
  obj_value = JuMP.objective_value(nonlinear_model)

  return status, x_sol, obj_value

end

```

Solve the problem to local optimality using `KNITRO`.

```julia
d = 10
status_KNITRO, x_sol_KNITRO, obj_value_KNITRO = solve_bilinear_locally_KNITRO(d; start_value_provided = :off, x_start_value = zeros(d), warm_start_type = :none, silent_flag = :off)
```
The output is:

```julia
=======================================
             Trial License
       (NOT FOR COMMERCIAL USE)
         Artelys Knitro 12.3.0
=======================================

Knitro presolve eliminated 0 variables and 0 constraints.

convex:                  0
datacheck:               0
hessian_no_f:            1
The problem is identified as a QCQP.

Problem Characteristics                                 (   Presolved)
-----------------------
Objective goal:  Maximize
Objective type:  quadratic
Number of variables:                                 10 (          10)
    bounded below only:                              10 (          10)
    bounded above only:                               0 (           0)
    bounded below and above:                          0 (           0)
    fixed:                                            0 (           0)
    free:                                             0 (           0)
Number of constraints:                               11 (          11)
    linear equalities:                                0 (           0)
    quadratic equalities:                             1 (           1)
    gen. nonlinear equalities:                        0 (           0)
    linear one-sided inequalities:                    1 (           1)
    quadratic one-sided inequalities:                 9 (           9)
    gen. nonlinear one-sided inequalities:            0 (           0)
    linear two-sided inequalities:                    0 (           0)
    quadratic two-sided inequalities:                 0 (           0)
    gen. nonlinear two-sided inequalities:            0 (           0)
Number of nonzeros in Jacobian:                      38 (          38)
Number of nonzeros in Hessian:                       11 (          11)

Knitro using the Interior-Point/Barrier Direct algorithm.

  Iter      Objective      FeasError   OptError    ||Step||    CGits
--------  --------------  ----------  ----------  ----------  -------
       0    0.000000e+00   1.000e+00
      10    9.413565e-01   8.790e-02   8.971e-01   1.222e-01        0
      20    2.397843e+01   1.159e-05   9.642e-03   1.307e-02        0
      25    2.400000e+01   6.485e-07   9.957e-08   1.186e-03        0

EXIT: Locally optimal solution found.

Final Statistics
----------------
Final objective value               =   2.39999992955648e+01
Final feasibility error (abs / rel) =   6.48e-07 / 6.48e-07
Final optimality error  (abs / rel) =   9.96e-08 / 1.99e-08
# of iterations                     =         25
# of CG iterations                  =          9
# of function evaluations           =          0
# of gradient evaluations           =          0
# of Hessian evaluations            =          0
Total program time (secs)           =       0.029 (     0.000 CPU time)
Time spent in evaluations (secs)    =       0.000
```

Note that `KNITRO` is providing the same objective value as `Gurobi`, but is not able certify the global optimality.

### Some helpful online links on speeding up the solution process

1. [How do I use MIP starts](https://support.gurobi.com/hc/en-us/articles/360043834831-How-do-I-use-MIP-starts-)
2. [Warmstart or VarHintVal in QP Python API](https://support.gurobi.com/hc/en-us/community/posts/360075269952-Warmstart-or-VarHintVal-in-QP-Python-API)
3. [Most important parameters for MIP models](https://www.gurobi.com/documentation/9.1/refman/mip_models.html)
