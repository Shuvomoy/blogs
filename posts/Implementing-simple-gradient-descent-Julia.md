@def title = "Simple gradient descent algorithm in Julia"
@def published = "April 10, 2020"
@def tags = ["programming", "Julia"]

# Implementing a simple gradient descent algorithm in Julia

**Shuvomoy Das Gupta**

*April 10, 2020*

In this blog, we discuss how to implement a simple gradient descent scheme in `Julia`. To do this, we will use the Julia package `ProximalOperators`, which is an excellent package to compute proximal operators and gradient of common convex functions. I highly recommend the package for anyone interested in operator splitting algorithms. You can find more information about the package at: [https://github.com/kul-forbes/ProximalOperators.jl](https://github.com/kul-forbes/ProximalOperators.jl).

---

**Table of contents**

\toc

---

**Jupyter notebook for this blog.** The jupyter notebook for this blog can be downloaded from [this link](https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/codes/implementing_simple_gradient_descent_Julia.ipynb) and viewed [here](https://github.com/Shuvomoy/blog/blob/gh-pages/codes/implementing_simple_gradient_descent_Julia.ipynb).  

Before we implement gradient descent method, we first record some necessary background.

### Background.

Given a differentiable convex function $f$, our goal is to solve $\textrm{minimize}\;f(x)$, where $x\in \mathbf{R}^n$ is the decision variable. To solve the problem above, we consider gradient descent algorithm. The gradient descent implements the following iteration scheme:

$
x_{n+1}  =  x_{n}-\gamma_{n}{\nabla f(x_{n})},\qquad (1)
$

where ${\nabla f(x_{n})}$ denotes a gradient of $f$ evaluated at the iterate $x_{n}$, and $n$ is our iteration counter. As our step size rule, we pick a sequence that is square-summable but not summable, e.g., $\gamma_{n}=1/n$, will do the job. 

We will go through the following steps:
1. Load the packages
2. Create the types
3. Write the functions

### Load the packages

Let us load the necessary packages that we are going to use.


```julia
## Load the packages to be used
# -----------------------------
# comment the first two lines if you already have ProximalOperators
using Pkg
Pkg.add("ProximalOperators")
using ProximalOperators, LinearAlgebra
```

###  Create the types

Next, we define a few Julia types, that we require to write an optimization solver in `Julia`. 

###  `GD_problem`

This type contains information about the problem instance, this bascially tells us what function $f$ we are trying to optimize over, one initial point $x_0$, and what should be the beginning step size $\gamma_0$.


```julia
struct GD_problem{F <: ProximableFunction,A <: AbstractVecOrMat{<:Real}, R <: Real}
    
    # problem structure, contains information regarding the problem
    
    f::F # the objective function
    x0::A # the intial condition
    γ::R # the stepsize
    
end
```

**Usage of `GD_problem`.**  For example, the user may wish to solve a simple least-squares problem using gradient descent. Then he can create a problem instance. A list of functions that we can use in this regard can be found in the documentation of `ProximalOperators`: [https://kul-forbes.github.io/ProximalOperators.jl/latest](https://kul-forbes.github.io/ProximalOperators.jl/latest).


```julia
# create a problem instance
# ------------------------

A = randn(6,5)

b = randn(6)

m, n = size(A)

# randomized intial point:

x0 = randn(n)

f = LeastSquares(A, b)

γ = 1.0

# create GD_problem

problem = GD_problem(f, x0, γ)
```



```julia 
    GD_problem(description : Least squares penalty
    domain      : n/a
    expression  : n/a
    parameters  : n/a, [1.0193204973421695, 1.168794122203917, 1.296856638205154, 0.5001687528217552, 0.1500452441023732], 1.0)
```



###  `GD_setting`

This type contains different parameters required to implement our algorithm, such as, 

* the initial step size $\gamma$, 
* maximum number of iterations $\textrm{maxit}$, 
* what should be the tolerance $\textrm{tol}$ (i.e., if $\| \nabla{f(x)} \| \leq \textrm{tol}$, we take that $x$ to be an optimal solution and terminate our algorithm), 
* whether to print out information about  the iterates or not controlled by a boolean variable $\textrm{verbose}$, and 
* how frequently to print out such information controlled by the variable $\textrm{freq}$.

The user may specify what values for these parameters above should be used. But if he does not specify anything, we should be able to have a default set of values to be used. We can achieve this by creating a simple constructor function for `GD_setting`.


```julia
struct GD_setting
    
    # user settings to solve the problem using Gradient Descent
    
    γ::Float64 # the step size
    maxit::Int64 # maximum number of iteration
    tol::Float64 # tolerance, i.e., if ||∇f(x)|| ≤ tol, we take x to be an optimal solution
    verbose::Boolean # whether to print information about the iterates
    freq::Int64 # how often print information about the iterates

    # constructor for the structure, so if user does not specify any particular values, 
    # then we create a GD_setting object with default values
    function GD_setting(; γ = 1, maxit = 1000, tol = 1e-8, verbose = false, freq = 10)
        new(γ, maxit, tol, verbose, freq)
    end
    
end
```

**Usage of `GD_setting`.** For the previously described least squares problem, we create the following setting instance.


```julia
setting = GD_setting(verbose = true, tol = 1e-2, maxit = 1000, freq = 100)
```

```julia 
    GD_setting(1, 1000, 0.01, true, 100)
```

###  `GD_state`

Now we define the type named  `GD_state` that describes the state our algorithm at iteration number $n$. The state is controlled by 

* current iterte $x_n$,
* the gradient of $f$ at the current iterate: ${\nabla{f}(x_n)}$,
* the stepsize at iteration $n$: $\gamma_n$, and
* iteration number: $n$.


```julia
mutable struct GD_state{T <: AbstractVecOrMat{<: Real}, I <: Integer, R <: Real} # contains information regarding one iterattion sequence
    
    x::T # iterate x_n
    ∇f_x::T # one gradient ∇f(x_n)
    γ::R # stepsize
    n::I # iteration counter
    
end
```

Also, once the user has given the problem information by creating a problem instance ``GD_problem``, we need a method to construct the initial value of the type `GD_state`,  as we did earlier for the least-squares problem. We create the initial state from the problem instance by writing a constructor function.


```julia
function GD_state(problem::GD_problem)
    
    # a constructor for the struct GD_state, it will take the problem data and create one state containing all 
    # the iterate information, current state of the gradient etc so that we can start our gradient descent scheme
    
    # unpack information from iter which is GD_iterable type
    x0 = copy(problem.x0) # to be safe
    f = problem.f
    γ = problem.γ
    ∇f_x, f_x = gradient(f, x0)
    n = 1
    
    return GD_state(x0, ∇f_x, γ, n)
    
end
```

```julia 
GD_state
```

## Write the functions

Now that we are done defining the types, we can now focus on writing the functions that will implement our gradient descent scheme. 

### `GD_iteration!`

First, we need a function that will take the problem information and the state of our algorithm at iteration number $n$, and then compute the next state for iteration number $n+1$ according to (1). 


```julia
function GD_iteration!(problem::GD_problem, state::GD_state)
    
    # this is the main iteration function, that takes the problem information, and the previous state, 
    # and create the new state using Gradient Descent algorithm
    
    # unpack the current state information
    x_n = state.x
    ∇f_x_n = state.∇f_x
    γ_n = state.γ
    n = state.n
    
    # compute the next state
    x_n_plus_1 = x_n - γ_n*∇f_x_n
    
    # now load the computed values in the state
    state.x = x_n_plus_1
    state.∇f_x, f_x = gradient(problem.f, x_n_plus_1) # note that f_x is not used anywhere
	# gradient(f,x) is a function in the ProximalOperators package, see its documentation 
	# if more information is required
    state.γ = 1/(n+1)
    state.n = n+1
    
    # done computing return the new state
    return state
    
end
```


```julia 
GD_iteration! (generic function with 1 method)
```

### `GD_solver`

Now we are in a position to write the main solver function named `GD_solver` that will be used by the end user. Internally, this function will take the problem information and the problem setting, and then it will

* create the initial state,
* keep updating the state using `GD_iteration!` function until we reach the termination criterion or the maximum number of iterations,
* print state of the algorithm if `verbose` is `true` at the specified frequency, and 
* return the final state.


```julia
## The solver function

function GD_solver(problem::GD_problem, setting::GD_setting)
    
    # this is the function that the end user will use to solve a particular problem, internally it is using the previously defined types and functions to run Gradient Descent Scheme
    # create the intial state
    state = GD_state(problem::GD_problem)
    
    ## time to run the loop
    while  (state.n < setting.maxit) & (norm(state.∇f_x, Inf) > setting.tol)
        # compute a new state
        state =  GD_iteration!(problem, state)
        # print information if verbose = true
        if setting.verbose == true
            if mod(state.n, setting.freq) == 0
                @info "iteration = $(state.n) | obj val = $(problem.f(state.x)) | gradient norm = $(norm(state.∇f_x, Inf))"
            end
        end
    end
    
    # print information regarding the final state
    
    @info "final iteration = $(state.n) | final obj val = $(problem.f(state.x)) | final gradient norm = $(norm(state.∇f_x, Inf))"
    return state
    
end
```


```julia 
GD_solver (generic function with 1 method)
```

​    



**Usage of `GD_solver`.** For the previously created `problem` and `setting`, we run our `GD_solver` function as follows.



```julia
# The following function will run the entire loop over the struct GradientDescent
```


```julia
final_state_GD = GD_solver(problem, setting)
```

```julia 
    ┌ Info: iteration = 100 |  
    │                 obj val = 0.10616975436389622 | 
    │                 gradient norm = 0.02964822736954073
    └ @ Main In[9]:16
    ┌ Info: iteration = 200 |  
    │                 obj val = 0.10526904472309405 | 
    │                 gradient norm = 0.019909231214547213
    └ @ Main In[9]:16
    ┌ Info: iteration = 300 |  
    │                 obj val = 0.10499388213419013 | 
    │                 gradient norm = 0.01578344670011611
    └ @ Main In[9]:16
    ┌ Info: iteration = 400 |  
    │                 obj val = 0.1048632843417343 | 
    │                 gradient norm = 0.013388371065530952
    └ @ Main In[9]:16
    ┌ Info: iteration = 500 |  
    │                 obj val = 0.10478781849295837 | 
    │                 gradient norm = 0.011784780173541287
    └ @ Main In[9]:16
    ┌ Info: iteration = 600 |  
    │                 obj val = 0.10473897267694926 | 
    │                 gradient norm = 0.010618620536011203
    └ @ Main In[9]:16
    ┌ Info: final iteration = 667 | 
    │     final obj val = 0.10471495099365061 | 
    │     final gradient norm = 0.009995364642054971
    └ @ Main In[9]:25





    GD_state([-0.747429923239562, 0.7857409304268766, -1.9752613372498569, -0.8109189160513681, 0.6537705422543739], [-0.002114026011632561, -0.008995100759876307, 0.009995364642054971, 0.0024806200462658134, -0.004665182937515999], 0.0014992503748125937, 667)
```

```julia
println("objective value found by our gradient descent $(f(final_state_GD.x))")

println("real objective value $(f(pinv(A)*b)) ")
```

```julia 
    objective value found by our gradient descent 0.10471495099365061
    real objective value 0.10452813292628083 
```

So, we do decent in terms of finding a good solution!

