@def title = "How to create a simple package in Julia"
@def published ="March 22, 2020"
@def tags =["programming", "Julia"]


# How to create a simple package in Julia
**Shuvomoy Das Gupta**

*March 22, 2020*

Thomson's Rule for First-Time Telescope Makers: *"It is faster to make a four-inch mirror then a six-inch mirror than to make a six-inch mirror."*

In this blog, we discuss how to create a simple package in `Julia` that computes the proximal operator of a convex quadratic function. The goal is just to illustrate how to create a structured solver package, the code underlying the solver is very inefficient. Let us call our package `TestPackage`. 

---

**Table of contents**

\toc

---

We will go through the following 6 steps in order to create a fully functioning package in Julia:

* Step 0. The inner mechanism of the package
* Step 1. Ensure `Git` and `Github` are set up properly
* Step 2. Install a few helpful packages for package development
* Step 3. Create the blank package `TestPackage` using `PkgTemplates`
* Step 5. Test the package and upload the files to Github

Now let us take a look at the steps in detail.

### Step 0. The inner mechanism of the package

In this step, we take a look at what our package is intended to do. We are given a quadratic function
$f(x) = \frac{1}{2}x^{T}Ax+b^{T}x+c$, where $A \in \mathbf{R}^{n \times n}$ is a positive semidefinite matrix, $ b\in\mathbf{R}^{n} $ is a column vector, and $ c\in\mathbf{R} $ is a scalar. Our goal is to compute the proximal map of the function $f$ at a point $v\in\mathbf{R}^{n}$ for the parameter $\gamma\in\mathbf{R}_{++}$, which is given by the following formula:

\begin{equation}
\mathbf{prox}_{\gamma f}(v)=(I+\gamma A)^{-1}(v-\gamma b).
\end{equation}

#### Structure of the package

We can structure our package in two parts: source code in the `src` folder, and a bunch of test files in the `test` folder. We will use a julia package named `PkgTemplates` that will create the aforementioned structure automatically for us. For now, let us just briefly take a look at the skeleton of our code.

**1. `src` folder**: The `src` folder has the following parts.

* *`Types.jl`* In this julia file, we store different data types to represent our problem data, parameter, and the output. The problem data essentially comprises of $A,b,c$. We will construct a data type `testPackageProblem` that will contain these matrices. Then, in the proximal map computation we have the parameter $\gamma$, which can be seen as a setting for our solver; we will create a data type called `testPackageSetting` that will contain this. Finally, our output vector will be stored in a data type named `testPackageResult`.

* *`Utils.jl`* This file contains a function for matrix inversion using the dependecy package `LinearAlgebra`.

* *`TestPackage.jl`* This is the main julia file corresponding to our package, indicated by its name. It will include both `Types.jl` and `Utils.jl`, and contain the main solver function called `solver_prox_quad`.

**2. `test` folder**: The test folder has file called `runtests.jl`, where we run our solver to different test datasets and see if everything is working as expected.

### Step 1. Ensure `Git` and `Github` are set up properly

Our package is going to reside in `Github` for version control purpose, so we need to create a free account first. Also, for version control, we need to have `git` installed on our computer as well. A very good instruction on how to perform these steps can be found [here](https://julia.quantecon.org/more_julia/version_control.html).

### Step 2. Install a few helpful packages for package development

Two very nice packages that can facilitate our package development by a significant margin are `PkgTemplates` and `Revise`. We install both packages first from Julia REPL.

```julia
using Pkg
```


```julia
Pkg.add("PkgTemplates")
```


```julia
Pkg.add("Revise")
```




A very convenient thing to do is to load `Revise` during startup. For this, please follow the instructions at this [link](https://timholy.github.io/Revise.jl/stable/config/#Using-Revise-by-default-1).

### Step 3. Create the blank package `TestPackage` using `PkgTemplates`

```julia
using PkgTemplates
```




Now create a specific template for our project.

```julia
template = Template(user = "shuvomoy") # change it to your specification
```




Now let us generate the package (empty for now).

```julia
generate("TestPackage", template)
```



As seen in the output above, the package is located at the folder `C:\Users\shuvo\.julia\dev\TestPackage`, copy this location; we are going to use it soon.



**Adding the project to the Julia Package Manager.** We need to go through the following steps, so that Julia's package manager knows about our project.



First, in the Julia REPL type the following: 

```julia

cd("C:\\Users\\shuvo\\.julia\\dev\\TestPackage")
# alternatively The following will also work
# cd(joinpath(DEPOT_PATH[1], "dev", "TestPackage"))
```

```julia
] activate # this will get us into the main Julia environment
```

```julia
] dev . # this adds our package TestPackage
```

```julia
] st # this is to see the change in our default package list
```



### Step 4. Create the julia files

#### Files in `src` folder ##

**Description of `Utils.jl`**

Create a file named `Utils.jl` in the `src` folder and copy the following code into it.

```julia
# create a file named Utils.jl and copy the following code
import LinearAlgebra # note that we are using another package LinearAlgebra in our file,
# so we need to add this package to the list of dependencies for our package

function mat_inv(A)
    return LinearAlgebra.pinv(A)
end
```




Note that in the code above, we have used the package `LinearAlgebra`. So, we need to add this package to the list of dependencies for our package. We can do that as follows.

```julia
] activate TestPackage
```


```julia
] add LinearAlgebra
```




To quit the active environment and return to the base again we can run the following.

```julia
] activate
```




**Description of `Types.jl`**

Create a file named `Types.jl` in the `src` folder and copy the following code into it.

```julia
# testPackageResult represents the result
struct testPackageResult
    x::Array{Float64}
end

# testPackageProblem contains the problem data

mutable struct testPackageProblem
    A::Array{Float64,2}
    b::Array{Float64}
    c::Float64
    # constructor
    function testPackageProblem(A,b,c)
        new(A,b,c)
    end
end

# testPackageSetting contains the problem setting, i.e., the parameter γ as set by the user

struct testPackageSetting
    γ::Float64
    # constructor
    function testPackageSetting(;γ = 1)
        new(γ)
    end
end
```




**Description of `TestPackage.jl`**

Create `TestPackage.jl` file in the `src` folder, and copy the following code in the file. Note that we are including both `Types.jl` and `Utils.jl` and exporting the name of contents, so that when we invoke `using TestPackage`, we can use those functions.

```julia
module TestPackage

include("./Types.jl")

export testPackageResult, testPackageProblem, testPackageSetting # so that the user can uses the types from the main module, see how to use this types in runtest.jl

include("./Utils.jl")

export mat_inv # so that the we can use the functions from the main module

function solver_prox_quad(v, problem::testPackageProblem, setting::testPackageSetting)
    # extract the problem data
    A = problem.A
    b = problem.b
    c = problem.c
    # extract the problem settings
    γ = setting.γ
    # compute the proximal parameter
    mat_prox_quad = LinearAlgebra.pinv(LinearAlgebra.I + γ*A)
    prox_quad_v = mat_prox_quad*(v-γ*b)
    return prox_quad_v
end

export solver_prox_quad

end # module
```




#### Files in `test` folder

Now that we have created the basic files for our package in the `src` folder, it is time now to test our code. To that goal, we create a file named `runtests.jl` in the `test` folder, where we have the following code.

```julia
#### Files in `test` folder ##
using TestPackage
using Test

@testset "testing for a certain set of data" begin
    # Write your own tests here.
    A = [ 2.67485     0.186124   0.583946  -0.145416   0.401269;
    -0.17813     0.526082  -0.325439  -1.46627   -3.09161 ;
    -1.52653    -0.631375  -0.618636  -2.36654    0.291429;
    -0.609349   -1.0533    -0.323403   1.34789   -0.409167;
    0.0687345  -0.946264  -0.118661  -1.49908   -1.30312]

    A = A*A'

    A = (A+A')/2

    b = [ 0.38052039048801956 ;
    0.5554436732038299 ;
    1.7334160203093987 ;
    0.4470254916964832 ;
    1.7676230141050555]

    c = 1.0

    v = ones(5)

    problem_data = testPackageProblem(A, b, c)

    problem_setting = testPackageSetting(γ = 2.0)

    result_approx = [0.10055689253917299, 0.24322623669390434, 0.0693976149433379, 0.09597235131516035, -0.530323239453701]

    @test solver_prox_quad(v, problem_data, problem_setting) ≈ result_approx

end

# we can potentially create more matrices like that
```




### Step 5. Test the package and upload the files to Github

Now we are in a position to test our package. We can do that in one of the following ways. We can just run the `runtests.jl` file from an editor such as Atom, and run it by pressing `shift-enter`.

Another way is open the Julia REPL, and run the following code.

```julia
] test TestPackage
```




So, everything is working fine!

**Add the package to Github.** Now, we just add our package to Git version control system. Go to your Github account, create a new repository with name `TestPackage.jl`. In the package creation page, please keep the boxes unchecked for the `README.md`, `LICENSE`, and `.gitignore`, because these are already taken care of when we created our package using `PkgTemplates`. Next, drag and drop the folder from our `~/.julia/dev` directory (to be specific, the contents of the folder `C:\Users\shuvo\.julia\dev\TestPackage`) to Github desktop. The final step is to click the `publish branch` button to upload our files to Github.

We are done creating a simple working Julia package!
