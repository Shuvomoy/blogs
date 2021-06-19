@def title = " Saving data files in Julia"
@def published ="April 28, 2020"
@def tags =["programming", "Julia"]


#  Saving data files in Julia
**Shuvomoy Das Gupta**

*April 28, 2020*


In this blog we will discuss ways of saving data into a file in Julia. We will discuss two methods of doing it: using `JLD2` package and writing directly into a `.jl` file. The literate code for this blog is available at: [this link](https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/codes/saving_data_files_julia.jl).

---

**Table of contents**

\toc

---

## Method 1: Using `JLD2` Package

Probably the quickest method to save the relevant data is to use the package `JLD2`.

```julia
using JLD2

A = randn(5,5)

b = randn(5)

str = "hello world"
```

Suppose, we want to save the variables `A`, `b`, and `str`, so that we could load them later for some other purpose.

Saving is done by using the following command.

```julia
@save "example.jld2" A b str
```

Now suppose we want to load the variables`A`, `b`, and `str` into our workspace. Then we can do that using the following command.

```julia
@load "example.jld2" A b str
```

Which will load the desired variables in our work space.

## Method 2: Writing into a `.jl` File
Now, in some circles, people prefer to save their data into some sort of text file that can be opened and seen using notepad etc. In such case, we could proceed as follows.

```julia
output_file = open("output_file.jl","w") # this will create a file named output_file.jl, where we will write the data.

write(output_file, "# We are saving the variables A, b, and str \n \n") # this line will act as a comment in the original file
```

The following command will write `A` into `output_file.jl`

```julia
write(output_file, "A = ") # writes A =

show(output_file, A) # writes the content of A

write(output_file, "; \n \n") # puts a semicolon to suppress the output and two line breaks
```

The following command will write `b` into `output_file.jl`

```julia
write(output_file, "b = ") # writes b =

show(output_file, b) # writes the content of b

write(output_file, "; \n \n") # puts a semicolon to suppress the output and two line breaks
```

The following command will write `str` into `output_file.jl`

```julia
write(output_file, "str = ") # writes str =

show(output_file, str) # writes the content of str

write(output_file, "; \n \n") # puts a semicolon to suppress the output and two line breaks
```

Finally we close the file.

```julia
close(output_file)
```

The benefit of this method is that we can open `output_file.jl` as a normal julia file and observe the contents of `A`, `b`, `str` which are preferred by some people. If we open the file, it will look something similar to the following.

```
# We are saving the variables A, b, and c

A = [0.2660919202329622 0.3836364273856932 -0.0843680697224458 -0.41802388053627293 -0.3642969656325985; -1.6326433489594656 0.8264017202126758 0.822230872741246 -2.131477875645448 -0.47154398683850046; 0.44945692740036236 1.8601083201122803 0.21328754042393291 1.6337810748330106 0.22377743588354782; 0.4650833154598731 0.7577147546636607 -0.3759919034593861 -0.7268603483224291 -1.355544453370908; 1.22775636169604 -0.6810273422844582 0.39118388831225326 0.6436653188719305 0.7072677389318383];

b = [-0.4480408521367389, -0.2470478314562805, -0.5170341188440505, -0.5941415056081019, -1.3602127485958821];

str = "hello world";
```

To load the variables in this case, we just run the following command.

```julia
include("output_file.jl")
```

**Some side note regarding generating the literate code:** First, cd to the directory where the file is in by running the a command similar to the following.
`cd("C:\\Users\\shuvo\\Desktop")`
Then, run
`using Literate`
Finally, run
`Literate.markdown("name_of_the_file.jl", "."; documenter=false)`

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

