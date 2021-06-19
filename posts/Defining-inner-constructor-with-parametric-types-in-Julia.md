@def title = " Defining inner constructor with parametric types in Julia"
@def published ="January 12, 2021"
@def tags =["programming", "Julia"]

#  Defining inner constructor with parametric types in Julia
**Shuvomoy Das Gupta**

*January 12, 2021*

Suppose we define the following parametric type in Julia. 

```julia 
abstract type CombinatorialProblem end

struct BinIntOptSF{V <: AbstractVector{<: Real}, M <: AbstractMatrix{<: Real}} <: CombinatorialProblem # means the concrete type BinIntOptSF is a subtype of the abstract type CombinatorialProblem
    c::V
    A::M
    b::V
end
```

This is a parametric type. We can create an instance of this type as follows. 

```julia 
m = 4
n = 5
c = randn(n)
A = randn(m,n)
b = randn(m)
sf_bin_int_opt = BinIntOptSF(c, A, b)
```

Now suppose, we want to enforce the following check in an instance of the type: `A` must have number of rows less than or equal to number of columns, number of rows of `A` is equal number of entries in `b`, and number of columns of `A` is equal to the number of entries in `c`. 

In this case we can enforce the check by writing the following inner constructor for the parametric type. 

```julia 
struct BinIntOptSF2{V <: AbstractVector{<: Real}, M <: AbstractMatrix{<: Real}} <: CombinatorialProblem

    c::V
    A::M
    b::V
    
    # Inner constructor 
    # -----------------
    function BinIntOptSF2{V,M}(c::V, A::M, b::V) where {V <: AbstractVector{<: Real}, M <: AbstractMatrix{<: Real}}
        
        if size(A,1) != size(b,1)
            error("A and b have incompatible direction")
        end
        
        if size(A,1) > size(A,2)
            error("A must have number of row less than or equal to number of columns")
        end
        
        if size(A,2) != size(c,1)
            error("A and c have incompatible dimension")
        end
        
        new{V, M}(c,A,b)
    end

# Writing the inner constructor above is going to replace the default constructor that is automatically defined by Julia. This may cause slight inconvenience, as we have to instantiate by running: bin_opt_instance = BinIntOptSF2{.,.}(...). Fortunately we can avoid this by adding the following line. 
    
    # Ensuring the one shot instantiation
    # -----------------------------------
    BinIntOptSF2(c::V, A::M, b::V) where {V <: AbstractVector{<: Real}, M <: AbstractMatrix{<: Real}} = BinIntOptSF2{V,M}(c,A,b) 

# We need to add this line so that we can instantiate by the convenient: bin_opt_instance = BinIntOptSF2(.). If we do not add this line, but call BinIntOptSF2(.), then we will get an error "No method matching...".
    
end
```

Let us instantiate the type above.

```julia 
sf_bin_int_opt_2 = BinIntOptSF2(c, A, b)
```

