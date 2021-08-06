@def title = "Running Julia code on MIT Supercloud"
@def published ="January 24, 2020"
@def tags =["programming", "Julia", "optimization"]

# Running Julia code on MIT Supercloud

**Shuvomoy Das Gupta**

*January 24, 2020*

In this blog, we are going to discuss how to MIT supercloud to run simple Julia code. 

---

\toc

---

As an illustrative example, I will use Julia code that uses lasso to approximately solve a sparse regression problem. 

First, we start with some common code that comes handy while working with the supercloud.

## Entering the supercloud

Assuming that we have setup our supercloud account already by following the instruction [here](https://supercloud.mit.edu/requesting-account), enter the supercloud from `bash` using the code:

`ssh tim@txe1-login.mit.edu
`
where replace `tim` with your username. 

## Transferring files from local computer to supercloud

If we want to transfer all the files in a folder called `PC_tests` to a folder named `supercloud_tests` located in the supercloud, we can achieve that by running the following command in `bash`.

`scp -r TIMs_Local_Folder_path/PC_tests/* tim@txe1-login.mit.edu:/home/gridsan/tim/supercloud_tests`

## Julia code

The code for the julia file is given below, please save it in a text file and name it ``lasso.jl``.


```julia
## Please copy the code in this code block and save it as lasso.jl
println("lasso.jl is running")
println("********************")

## The following function is used to compute cardinality of a vector.
function l_0_norm(x,y,n,delta_l0_norm)
    diffxy = x-y
    l_0_norm_xy = ones(n)
    for i in 1:n
        if abs(diffxy[i]) <= delta_l0_norm
            l_0_norm_xy[i]=0
        end
    end
    return sum(l_0_norm_xy)
end

function card(x, n, delta_l0_norm)
    return l_0_norm(x,zeros(n),n,delta_l0_norm)
end

## Data

b = randn(100)

A = randn(100,200)

m, n = size(A)

λ_glment_array = 10 .^(range(-4,stop=4,length=1000))

k = 20

## running lasso using glmnet
using GLMNet, LinearAlgebra

path = glmnet(A, b, lambda = λ_glment_array)

matrix_x = convert(Matrix, path.betas)

## finding glmnet smallest

λ_glmnet_smallest = 0

x_glmnet = zeros(n)

m1, n1 = size(matrix_x)

for i in 1:n1
    cardinality_present = card(matrix_x[:,i], m1, 1e-5)
    println("cardinality of the ", i , "-th solution is = ", cardinality_present)
    if cardinality_present <= k
        println("card(signal) =", cardinality_present, " found for λ = ",  λ_glment_array[i])
        global λ_glmnet_smallest =  λ_glment_array[i]
        println(" corresponding solution x is = ", matrix_x[:,i])
        global x_glmnet = matrix_x[:,i]
        break
    end
end

λ_smallest = 2*m*λ_glmnet_smallest

## check objecgive value
obj_glmnet = (norm(A*x_glmnet - b))^2

# finding the index of nonzero elements
Z = findall(iszero, x_glmnet)
k == n - length(Z)

##
using SCS
solver = SCSSolver()
using Convex

x = Variable(n)
objective = sumsquares(A * x - b)
problem = minimize(objective)
for j = 1:length(Z)
    problem.constraints += [x[Z[j]] == 0.0]
end
solve!(problem, solver)
x_convex = x.value

function polish(x)
    n = length(x)
    x_polished = x
    for i in 1:n
        if abs(x[i]) <= 1e-5
            x_polished[i] = 0
        end
    end
    return x_polished
end

x_convex_polished = polish(x_convex)

obj_convex = norm(A*x_convex_polished - b)^2

using JLD2

println("saving the data now...")

@save "lasso_data.jld2" A b λ_glment_array x_glmnet x_convex_polished obj_glmnet obj_convex

println("done running the code :)")

```

## Shell script to submit the job

Now we are going to create a shell script that will be used to submit the job. The code for the shell script is below. Please save it in a text file, and name it ``run_lasso.sh``.


```julia
#!/bin/bash

#SBATCH -o lasso_results.txt
#SBATCH --nodes=1 # Number of node
#SBATCH --ntasks=1 # Number of tasks
#SBATCH --cpus-per-task=32 # How many threads to assign
#SBATCH --mem=32G # Hom much memory 

# Initialize the module command first source
source /etc/profile

# Load Julia Module
# Find out the list of modules available by running the command
#  module avail
module load julia/1.5.2

# Call your script as you would from the command line
julia lasso.jl

```

## Submitting the job

Now log in to MIT supercloud, copy the files created above to your working directory, and run the following command.


```julia
LLsub run_lasso.sh
```

## If we need more memory

Here we should keep in mind that MIT supercloud nodes have 16 cores and 64 GB of RAM in total. As a result, each core has about 4 GB. So if our job requires a certain amount of memory, we need to allocate memory accordingly. For example, if our program needs roughly 32 GB memory, then we need to allocate 8 cores, which can be done by running the following command.


```julia
LLsub myScript.sh -s 8
```

## Observing the status of the submitted job

After submitting the job, we can check the status of our jobs by running the following command


```julia
LLstat
```

which will have an output like:


```julia
LLGrid: txe1 (running slurm 19.05.5)
JOBID     ARRAY_J    NAME         USER     START_TIME          PARTITION  CPUS  FEATURES  MIN_MEMORY  ST  NODELIST(REASON)
17019     40986      run_lasso       Student  2020-10-19T15:35:46 normal     1     xeon-e5   5G          R   gpu-2
17214     40980      myJob1       Student  2020-10-19T15:35:37 normal     1     xeon-e5   5G          R   gpu-2

```

## Killing a job

If we want to kill one of the jobs, e.g., ``run_lasso``, then we can do that by running the following command:


```
LLkill 17019
```

## Installing updated Julia package

Sometimes due to conflicts, one may have trouble installing an updated version of a Julia package.  In that case run the following code. 

First, on bash:  (if necessary, you can add the following to the `.bashrc` file)

```bash
export TMPDIR=/home/gridsan/tim/TMPDIR/ 
# change the default /tmp directory, ensure that the   
# TMPDIR folder exists in the user folder exists
# you can check if the folder changed by running 
# julia>  tempdir()

unset JULIA_LOAD_PATH 
# remove the shared version of julia load path
# you check the JULIA_LOAD_PATH by running 
# julia> LOAD_PATH

export JULIA_DEPOT_PATH=/home/gridsan/tim/.julia/ 
# remove the shared 
# version of julia depot path
# you check the JULIA_LOAD_PATH by running 
# julia> DEPOT_PATH

julia # start julia
```

Then in Julia run:

```julia 
 add JLD2@0.4.11 # or any version that you want
```

Going forward one need not run these extra steps, just loading the Julia module suffices.

## Useful link

A very good link for understanding slurm scripts (how to write the .sh file) is:

[https://researchcomputing.princeton.edu/support/knowledge-base/slurm](https://researchcomputing.princeton.edu/support/knowledge-base/slurm)

