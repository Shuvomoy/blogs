@def title = "Tips and tricks in Mathematica"
@def published ="August 26, 2021"
@def tags =["programming", "Mathematica"]

# Tips and tricks in Mathematica

**Shuvomoy Das Gupta**

*August 26, 2021*

In this blog, I am collecting a bunch of tips and tricks in Mathematica, that I find very useful for my work.

---

\toc 

---

## Clear all data

```mathematica
(*Important simplification codes*)
ClearAll["Global`*"];
```



## Removing outliers from data

This solution is due to [Carl Lange](https://mathematica.stackexchange.com/users/57593/carl-lange).

```mathematica
(*This dataset contains some outliers*)
Data={{0.105, 0.989213}, {0.106414, 0.988926}, {0.107828, 
  0.988636}, {0.109242, 0.988343}, {0.110657, 0.988049}, {0.112071, 
  0.987748}, {0.113485, 0.}, {0.114899, 1.}, {0.116313, 
  0.986826}, {0.117727, 0.986512}, {0.119141, 0.986196}, {0.120556, 
  0.995073}, {0.12197, 0.985551}, {0.123384, 0.0154883}, {0.124798, 
  0.984894}, {0.126212, 1.}, {0.127626, 0.984222}, {0.12904, 
  0.983887}, {0.130455, 0.983538}, {0.131869, 0.983197}, {0.133283, 
  0.}, {0.134697, 0.970927}, {0.136111, 0.98213}, {0.137525, 
  0.98177}, {0.138939, 1.}, {0.140354, 0.981041}, {0.141768, 
  0.980672}, {0.143182, 0.826229}, {0.144596, 0.979923}, {0.14601, 
  0.979546}, {0.147424, 0.979163}, {0.148838, 0.978778}, {0.150253, 
  0.978392}, {0.151667, 0.978}, {0.153081, 0.977605}, {0.154495, 
  0.977208}, {0.155909, 0.976807}, {0.157323, 0.976404}, {0.158737, 
  0.975999}, {0.160152, 0.55766}, {0.161566, 
  0.975177}, {0.16298, -0.000401533}, {0.164394, 0.974344}, {0.165808,
   1.00182}, {0.167222, 0.}, {0.168636, 0.973073}, {0.170051, 
  0.972646}, {0.171465, 0.972211}, {0.172879, 0.971787}, {0.174293, 
  0.971338}, {0.175707, 0.970898}, {0.177121, 0.970455}, {0.178535, 
  0.97001}, {0.179949, 0.96956}, {0.181364, -0.000767749}, {0.182778, 
  0.968655}, {0.184192, 0.968197}, {0.185606, 0.967738}, {0.18702, 
  0.967275}, {0.188434, 0.96681}, {0.189848, 0.966343}, {0.191263, 
  0.}, {0.192677, 0.965404}, {0.194091, 0.964925}, {0.195505, 
  0.964447}, {0.196919, 0.963967}, {0.198333, 0.963484}, {0.199747, 
  0.962999}, {0.201162, 1.}, {0.202576, 0.962022}, {0.20399, 
  0.961529}, {0.205404, 0.961034}, {0.206818, 0.960536}, {0.208232, 
  0.960036}, {0.209646, 0.959534}, {0.211061, 0.959029}, {0.212475, 
  0.958522}, {0.213889, 0.958013}, {0.215303, 1.}, {0.216717, 
  0.956987}, {0.218131, 0.956471}, {0.219545, 0.955953}, {0.22096, 
  0.955432}, {0.222374, 0.954909}, {0.223788, 0.954385}, {0.225202, 
  0.894605}, {0.226616, 0.953327}, {0.22803, 0.952796}, {0.229444, 
  0.952262}, {0.230859, 0.951726}, {0.232273, 0.951188}, {0.233687, 
  0.950648}, {0.235101, 0.950106}, {0.236515, 0.949561}, {0.237929, 
  0.949017}, {0.239343, 0.948467}, {0.240758, 0.947917}, {0.242172, 
  0.947364}, {0.243586, 0.946811}, {0.245, 0.946254}};
  
(*Remove outliers*)
CleanedData = DeleteAnomalies[
  LearnDistribution[MovingMedian[Data, 5], Method -> "Multinormal"], 
  Data]

(*Plot the cleaned data*)
ListPlot[CleanedData, PlotRange -> Full]

(*Compare with original data, red is original data, blue is cleaned data*)
ListPlot[{Data, CleanedData}, PlotStyle -> {Red, Blue}]
```

![image-20210824081939805](https://raw.githubusercontent.com/Shuvomoy/blogs/master/posts/Tips_and_tricks_in_Mathematica.assets/image-20210824081939805.png)

![image-20210824081952993](https://raw.githubusercontent.com/Shuvomoy/blogs/master/posts/Tips_and_tricks_in_Mathematica.assets/image-20210824081952993.png)

## Collecting terms with specific pattern

```mathematica
(*This code will collect terms with a specific patterns*)
(*Caution all the terms have to be scalrs, does not work with 
table term such x[i] etc, but works with xi and so on*)

CollectWRTVarList[expr_, vars_List] := 
  Expand[Simplify[
     expr /. Flatten[
       Solve[# == ToString@#, First@Variables@#] & /@ vars]], 
    Alternatives @@ ToString /@ vars] /. 
   Thread[ToString /@ vars -> vars];
            
(*Example*)
            
CollectWRTVarList[
 a c x1 + a d x1 + a c x2 + a d x2 + a y1 + a y2, {x1 + x2, y1 + y2}]
            
(*output = a (c + d) (x1 + x2) + a (y1 + y2)*)            
```



## Finding simple seqeunce function from the data

```mathematica
(*Example 1.*)
d[0] = d0;
d[i_] := 2 L + d[i - 1] /; i >= 1 

(*here/;i\[GreaterEqual]1 imposes the condition that i 
should be greater than 1 for the functional description to be true*);

dAnalytic[n_] := FindSequenceFunction[Table[d[i], {i, 1, 10}], n]

(*Example 2.*)
b[0] = b0;
b[i_] := 2 + (d[i - 1]/L) + b[i - 1] /; i >= 1
bAnalytic[n_] := FindSequenceFunction[Table[b[i], {i, 1, 10}], n]

(*CAUTION:the output is 0-based formula if in the seqeuence function 
we have {i,1,10},i.e.,in dAnalytic[n] and bAnalytic[n],n will range 
from {0,1,...},if we had started {i,0,10} the formulas would have n 
ranging from {1,2,...}*)
                
StringForm[" d[n] = `1` and b[n] = `2` for n \[Element] {0,1,2...}", dAnalytic[n], bAnalytic[n]]                
```

The output is:

![image-20210824083117487](https://raw.githubusercontent.com/Shuvomoy/blogs/master/posts/Tips_and_tricks_in_Mathematica.assets/image-20210824083117487.png)

## Simplifying recursion without closed form solution

```mathematica
(*----------------------------------------------*)
w[0] = w0;
w[k_] := w[k - 1] - 1/L Sum[h[k, i] g[i], {i, 0, k - 1}] /; k >= 1 

(* here /; i\[GreaterEqual] 1 imposes the condition that i \
should be greater than 1 for the functional description to be true *);

StringForm["Table of w[i]"]
Table[w[i], {i, 0, 3}] // TableForm

(*Observe the pattern in w[i]*)

StringForm["Observe the pattern in w[i]"]
{Expand[w[1]],
 Collect[Expand[w[2]], {g[0], g[1]}],
 Collect[Expand[w[3]], {g[0], g[1], g[2]}]}
 
(*Try to guess the pattern now*)
StringForm["Guess the pattern in w[i] via wSimp[i]"]
wSimp[0] = w0;
wSimp[k_] := wSimp[0] - 1/L Sum[Sum[h[l, i] g[i], {l, i + 1, k}], {i, 0, k - 1}] /; k >= 1

Table[w[i], {i, 1, 3}] // TableForm
Table[Simplify[wSimp[i]], {i, 1, 3}] // TableForm

StringForm["See if w[i]==wSimp[i]"]
Table[Simplify[wSimp[i] - w[i]], {i, 1, 3}] // TableForm
```

![image-20210824084532271](https://raw.githubusercontent.com/Shuvomoy/blogs/master/posts/Tips_and_tricks_in_Mathematica.assets/image-20210824084532271.png)

## Define functions that will take vectors

```mathematica
(*Define functions that will take vectors*)
fun2[a_, b_] := 
 Assuming[a \[Element] Vectors[n, Reals] && 
    b \[Element] Vectors[n, Reals], (c1 a + c2 b) . (c1 a + c2 b) // 
     TensorExpand // Simplify] /. p_ Dot[t_, t_] -> p Norm[t]^2
fun2[x, y]

(*Test*)
fun2[x, y] 

(*Output = c1^2 Norm[x]^2 + c2 (2 c1 x . y + c2 Norm[y]^2) *)
```

## Running `Mathematica` from `Julia`

We can run Mathematica from Julia using `MathLink.jl`. We can install `MathLink.jl` in Julia using the following command

```julia 
ENV["JULIA_MATHKERNEL"]="C:\\Program Files\\Wolfram Research\\Mathematica\\12.3\\MathKernel.exe" # path to MathKernel.exe
ENV["JULIA_MATHLINK"]="C:\\Program Files\\Wolfram Research\\Mathematica\\12.3\\SystemFiles\\Links\\MathLink\\DeveloperKit\\Windows-x86-64\\SystemAdditions\\ml64i4.dll" # path to MathLink
import Pkg
Pkg.add("MathLink")
Pkg.build("MathLink")
```

Now suppose we want to run the following Mathematica command in Julia.

```mathematica
a = Table[{x, N[x Sin[x]]}, {x, 0, 4, .3}];
FindFormula[a, x]
```

We run this Julia by putting the expression above in `weval(W``)` block as follows.

```julia 
weval(W`
    a = Table[{x, N[x Sin[x]]}, {x, 0, 4, .3}];
    FindFormula[a, x]
    `)

# output
# ------
# W"Times"(W"x", W"Sin"(W"x"))
# In Mathematica it will correspond to
# Times[x,Sin[x]] = x Sin[x]
```

## Copy Mathematic code as Unicode characters

This solution is from the [link](https://mathematica.stackexchange.com/questions/1137/how-to-copy-as-unicode-from-a-notebook). Suppose we want to copy the following text from Mathematica into some other editor while preserving the Unicode characters:

```mathematica
(*Code to copy from Mathematica preserving Unicdoe*)
u1 = -((b + β + a1)/((b + β) rμ));
u2 = -((μ a2)/((1 + b μ + β μ) rμ)); 
Reduce[
 Abs[1 + rμ u1] > Abs[rμ u2] && 
  b > 0 && β > 0 && μ > 0 && 0 < rμ <= 1 && 
  Abs[u1] >= 1 && a1 > 0, {a1, a2}, Reals]
```

We define the following function.

```mathematica
SetAttributes[copyUnicode, HoldFirst]

copyUnicode[expr_, form_: InputForm] := 
  Run["clip <", 
   Export["$Clipboard.temp", ToString[Unevaluated@expr, form], "Text", 
    CharacterEncoding -> "Unicode"]];
```

And then run the following from Mathematica.

```mathematica
(*Put the entire Mathematica code in the function copyUnicode, i.e., (code_to_copy)//copyUnicode or copyUnicode[code_to_copy]*)
(u1 = -((b + β + a1)/((b + β) rμ));
u2 = -((μ a2)/((1 + b μ + β μ) rμ)); 
Reduce[
 Abs[1 + rμ u1] > Abs[rμ u2] && 
  b > 0 && β > 0 && μ > 0 && 0 < rμ <= 1 && 
  Abs[u1] >= 1 && a1 > 0, {a1, a2}, Reals])//copyUnicode    
  
(*output:
﻿Hold[u1 = -((b + β + a1)/((b + β)*rμ)); u2 = -((μ*a2)/((1 + b*μ + β*μ)*rμ)); 1*Reduce[Abs[1 + rμ*u1] > Abs[rμ*u2] && b > 0 && β > 0 && μ > 0 && Inequality[0, Less, rμ, LessEqual, 1] && Abs[u1] >= 1 && a1 > 0, {a1, a2}, Reals]]
*)  

```

## Parametric optimization problem

```julia 
Minimize[{-a x + (b + \[Beta])/2 (x^2 + y^2), 
  a > (b + \[Beta]) && \[Beta] > 0 && b > 0 && 0 <= x <= 1 && 
   0 <= y <= 1 }, {x, y}]
```



