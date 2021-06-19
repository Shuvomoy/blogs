@def title = "Mutating function and array assignment in Julia"
@def published ="January 18, 2021"
@def tags =["programming"]

# Mutating function and array assignment in Julia
**Shuvomoy Das Gupta**

*January 18, 2021*


A mutating function in Julia changes one of its inputs. A common convention is to use `!` after the name of the function if it is a mutating function. 

If we are mutating an array type object, for proper mutation we must change what the array contains, *not what the array points to*.

```julia
function frst_inpt_dbl_of_2nd_inpt!(A::AbstractMatrix{<: Real}, B::AbstractMatrix{<: Real})
    A .= 2*B # this .= notation ensures that we are changing the content of A
    # 💀 if we do A = B, then it will only change the reference to the array A
end
```
Let us test the function above.

```julia
B = ones(5,5)
A = similar(A)
frst_inpt_dbl_of_2nd_inpt!(A, B)
A
```
which indeed gives us the output:
```julia
 A =
 2.0  2.0  2.0  2.0  2.0
 2.0  2.0  2.0  2.0  2.0
 2.0  2.0  2.0  2.0  2.0
 2.0  2.0  2.0  2.0  2.0
 2.0  2.0  2.0  2.0  2.0
```
So the contents of `A` has indeed changed.

Just to experiment, let us experiment with a *wrong* implementation below.

```julia
function wrong_frst_inpt_dbl_of_2nd_inpt_💀!(A::AbstractMatrix{<: Real}, B::AbstractMatrix{<: Real})
    A = 2*B # this will change the references to A, not the values
end
```

```julia
B = ones(5,5)
A = similar(A)
wrong_frst_inpt_dbl_of_2nd_inpt_💀!(A, B)
println(A)
```

which gives us the disastrous output:
```julia
A =
3.5e-323  8.4e-323   1.43e-322  2.17e-322  2.8e-322
4.0e-323  9.4e-323   1.73e-322  2.4e-322   2.87e-322
4.4e-323  1.1e-322   1.9e-322   2.5e-322   3.1e-322
7.0e-323  1.14e-322  2.03e-322  2.57e-322  3.16e-322
8.0e-323  1.33e-322  2.08e-322  2.67e-322  3.2e-322
```

For a more detailed explanation of why this happens, see the blog by John White [here](http://www.johnmyleswhite.com/notebook/2014/09/06/values-vs-bindings-the-map-is-not-the-territory/).