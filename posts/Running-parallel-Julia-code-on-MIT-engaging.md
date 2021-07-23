@def title = "Running parallel Julia code on MIT Engaging"
@def published ="February 1, 2021"
@def tags =["programming", "Julia"]

# Running parallel Julia code on MIT Engaging
**Shuvomoy Das Gupta**

*February 1, 2021*

In this blog, we will discuss how to run  parallel `Julia` code on MIT Engaging. 

---

**Table of contents**

\toc

---

## Setup

We can achieve this task by using parallelization techniques provided in `Julia`. For this blog, we will consider `pmap`. In the following code, we create 10 large matrices and then perform a singular value decomposition on each. We will show that parallelizing this computation can attain significant speedup simply.

**Julia Code.** The code for the `Julia` file is given below. Please save it in a text file and name it `pmap_julia.jl`. 

```julia 
using ClusterManagers,Distributed, BenchmarkTools

# Add in the cores allocated by the scheduler as workers

addprocs(SlurmManager(parse(Int,ENV["SLURM_NTASKS"])-1))

print("Added workers: ")

println(nworkers())

using LinearAlgebra # Load the package centrally

@everywhere using LinearAlgebra # load it on each process

# DATA GENERATION: The following code will be run only on one core
# ----------------------------------------------------------------

# create the array of random initial points
X_array=[rand(100,100) for i in 1:10]

# benchmark serial implementation
b1 = @benchmark map(svd, X_array)

println("benchmark for serial code")
println("*************************")
io = IOBuffer()
show(io, "text/plain", b1)
s = String(take!(io))
println(s)

# benchmark parallel implementation
b2 = @benchmark pmap(svd, X_array)

println("benchmark for parallel code")
println("***************************")
io = IOBuffer()
show(io, "text/plain", b2)
s = String(take!(io))
println(s)
```

Note that in the first line, we have loaded the packages that we use. If they are not installed, you can install them by running the following commands from `bash`. 

```julia 
srun --pty -p sched_mit_sloan_interactive julia
```

This will start a `Julia REPL` in `bash`. Run the following commands: 

```julia 
] add ClusterManagers,Distributed, BenchmarkTools
```

## Shell script to submit the job

Now we are going to create a shell script that will be used to submit the job. The code for the shell script is below. Please save it in a text file, and name it `run_pmap_julia.sh`. In the code, `SBATCH -o pmap_julia.log-%j` indicates the name of the file where the output is written, and `SBATCH -n 14` indicates the number of cores or cpus allocated to the job. 

```julia 
#!/bin/bash

# Slurm sbatch options
#SBATCH -o pmap_julia_log-%j.txt
#SBATCH -n 14

# Initialize the module command first source
source /etc/profile

# Load Julia Module
module load julia/1.5.2
   
 
# Call your script as you would from the command line
# Call your script as you would from the command line
julia pmap_julia.jl
```

## Submitting the job

Now log in to MIT engaging, copy the files created above to your working directory, and run the following command. Ensure that the `pwd` command results in the working directory, else use the command `cd directory_that_contains_the_code` to change the working directory.

```
sbatch run_pmap_julia.sh
```

That's it! Once the computation is done, we see from the output log file that parallelization has decreased the computation time significantly: 

```julia 
benchmark for serial code
*************************
BenchmarkTools.Trial: 
  memory estimate:  4.71 MiB
  allocs estimate:  102
  --------------
  minimum time:     52.735 ms (0.00% GC)
  median time:      53.473 ms (0.00% GC)
  mean time:        53.792 ms (0.15% GC)
  maximum time:     56.662 ms (4.62% GC)
  --------------
  samples:          93
  evals/sample:     1
    
benchmark for parallel code
***************************
BenchmarkTools.Trial: 
  memory estimate:  1.60 MiB
  allocs estimate:  1436
  --------------
  minimum time:     6.624 ms (0.00% GC)
  median time:      6.927 ms (0.00% GC)
  mean time:        7.172 ms (0.78% GC)
  maximum time:     12.139 ms (0.00% GC)
  --------------
  samples:          697
  evals/sample:     1
```

where we see that the parallel code is more than 7 times faster than that of the serial code!

## Some handy commands

 The `sbatch` commands are the same as `MIT Supercloud`, so we can use the same command as mentioned in the blogs about the `MIT Sueprcloud`.

* To protect the data:

```julia 
eo-fix-storage
eo-fix-permissions all
```

* To see the storage: 

```julia 
eo-show-quota
```

* Note that, each user has access to 4 different storage areas:

1. /home/myusername - working space for source code, scripts,hand-edited files etc.  Each user has 100GB by default.
2. /pool001/myusername - Extra storage space.  Each user has 1 TB by default
3. /nobackup1/myusername - Very fast lustre parallel file system for parallel I/O, this should not be used for long-term storage.
4. /nfs/sloanlab001/projects - Requestable shared project space that can be linked to your home folder, eg $HOME/projects/myproject_proj

* To view the running jobs we can type the command:

```
eo-show-myjobs
```

* Suppose we want to cancel job number 12345. The command is:

```julia 
scancel 12345
```

* If we want to get a quick view of all the jobs completed within the last 5 days, we use:

```
eo-show-history 
```

## Helpful Link

A comprehensive documentation about `Engaging` is available at the link:

[https://wikis.mit.edu/confluence/display/sloanrc/Engaging+Platform](https://wikis.mit.edu/confluence/display/sloanrc/Engaging+Platform)

which requires an MIT login. 