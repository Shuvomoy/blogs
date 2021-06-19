@def title = "Verifying inequalities in Mathematica"
@def published ="November 15, 2020"
@def tags =["programming", "Mathematica"]


# Verifying inequalities in Mathematica
**Shuvomoy Das Gupta**

*November 15, 2020*

In applied mathematics, we often need to use inequalities to simplify our computation. Of course, verifying an inequality would requires picking up pen and paper and proving the it rigorously. However, a good idea prior to that proving phase is to test if the inequality holds for smaller dimensions. This verification for an inequality over smaller dimension can be done efficiently using `Mathematica`. Here are two simple examples.

\toc

### A simple example: AM-GM inequality

We start with a very simple example: the well-known AM-GM inequality. It states that for $x>0,y>0$ we have 
$$
\frac{x+y}{2} \geq \sqrt{xy}.
$$
We can verify it as follows. 

```mathematica
(* Clear all the variables, this often comes handy*)
ClearAll["Global`*"];

(*construct the conditions*)
conditions = x > 0 && y > 0;

(*check wheather the AM-GM inequality holds for all x,y satisfying conditions*)

(*Create the AM GM inequality*)
inequalityAMGM = 
 ForAll[{x, y}, 
  conditions, (x + y)/2 >= 
   Sqrt[x y]] 
   
(*Verify if the inequality holds for all x,y satisfying conditions*)   
Resolve[inequalityAMGM] 
```

where we get the output `True`, so we have verified the AM-GM inequality.

### Cauchy inequality

The Cauchy inequality probably one of the most famous inequalities. Let us verify it in `Mathematica` for dimension 3.

```mathematica 
(* Clear all the variables *)
ClearAll["Global`*"];

(* Create the Cauchy inequality *)
ineqCauchy = ForAll[{x, y}, Element[x | y, Vectors[3, Reals]], 
       Abs[x . y] <= Norm[x]*Norm[y]]; 

(* Verify if the inequality holds *)
Resolve[ineqCauchy]
```

which outputs `True` again. We can run this for larger dimension too. However, keep in mind that larger the dimension, longer it would take for `Mathematica` to verify it. Hence, it is best if this verification process is kept confined to a smaller dimension, and then if the verification process yields `True`, then go for the good old pen and paper to prove it formally. 

### A more complicated example: Performance Estimation Problem

As our final example, we consider verifying an inequality that shows up in the performance estimation problem, for more details about this problem, please see the paper by Taylor et al. [here](https://arxiv.org/pdf/1502.05666.pdf). We want to verify an identity that shows up in Theorem 4 of the mentioned paper. Given vectors $f_i,f_i,x_i,x_j,g_i,g_j\in \mathbf{R}^n$ we want to show that the following two terms $t_1$ and $t_2$ are equal to each other. We have:
$$
t_1 = f_i - f_j - \left\langle g_j \mid x_i - x_j \right\rangle - \frac{1}{2(1/(\mu / L))} \left( \frac{1}{L}\|g_i - g_j \|_2^2 + \mu \|x_i - x_j \|_2^2-\left\langle g_j - g_i \mid x_j - x_i\right\rangle\right)
$$


and


$$
\begin{align*}
s_{i} & =\frac{\mu}{L-\mu}\left\langle g_{i}\mid x_{i}\right\rangle -\frac{\mu L}{2(L-\mu)}\|x_{i}\|_{2}^{2}-\frac{1}{2(L-\mu)}\|g_{i}\|_{2}^{2}\\
s_{j} & =\frac{\mu}{L-\mu}\left\langle g_{j}\mid x_{j}\right\rangle -\frac{\mu L}{2(L-\mu)}\|x_{j}\|_{2}^{2}-\frac{1}{2(L-\mu)}\|g_{j}\|_{2}^{2}\\
p_{ij} & =\left\langle g_{j}-\mu x_{j}\mid\left(\frac{L}{L-\mu}x_{i}-\frac{1}{L-\mu}g_{i}\right)-\left(\frac{L}{L-\mu}x_{j}-\frac{1}{L-\mu}g_{j}\right)\right\rangle \\
t_{2} & =s_{i}-s_{j}-p_{ij}.
\end{align*}
$$
Where $0<\mu < L$.

First we clear the variables.

```mathematica
(*Clear all the variables*)
ClearAll["Global`*"];
```

We are going to do this test for $n=3$. Let us create our assumptions.

```mathematica
n = 3;

myAssumptions = (xi | xj | gi | gj) \[Element] 
   Vectors[n, Reals] && (\[Mu] | L) \[Element] Reals && \[Mu] > 0 && 
  L > 0 && \[Mu] < L
```

Let us construct $t_1$ first.

```mathematica
t1 = fi - fj - 
 gj.(xi - xj) - (-((2 \[Mu] (-gi + gj).(-xi + xj))/L) + 
  Norm[gi - gj]^2/L + \[Mu] Norm[xi - xj]^2)/(2 (1 - \[Mu]/L))
```

Now, let us construct $t_2$ by constructing $s_i,s_j,$ and $p_{ij}$.

```mathematica
si = \[Mu]/(L - \[Mu]) Dot[gi, 
     xi] - (\[Mu] L)/(2 (L - \[Mu])) Norm[xi]^2 - 
   1/(2 (L - \[Mu])) Norm[gi]^2 + fi;

sj = \[Mu]/(L - \[Mu]) Dot[gj, 
     xj] - (\[Mu] L)/(2 (L - \[Mu])) Norm[xj]^2 - 
   1/(2 (L - \[Mu])) Norm[gj]^2 + fj;
     
pij = Dot[
   gj - \[Mu] xj, (L/(L - \[Mu]) xi - 
      1/(L - \[Mu]) gi) - (L/(L - \[Mu]) xj - 1/(L - \[Mu]) gj)];
      
t2 = si - sj - pij
```

Construct the identity now. 

```mathematica
identityPEP = ForAll[{gi, gj, xi, xj, \[Mu], L}, myAssumptions, t1== t2]
```

Time to resolve.

```julia 
Resolve[identityPEP]
```

This yields `True` (it will take a few minutes) so the identity in question is true for $n=3$. 

