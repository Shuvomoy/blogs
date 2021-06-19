@def title = "Computing scaled relative graph of composition of operators in Mathematica"   
@def published = "November 17, 2020"   
@def tags =["programming", "Mathematica"]  


# Computing scaled relative graph of composition of operators in Mathematica

**Shuvomoy Das Gupta**

*November 17, 2020*


This blog is based on a question I posted [here](https://mathematica.stackexchange.com/questions/233264/plotting-minkowski-product-of-two-sets-in-complex-2d-plane) on [https://mathematica.stackexchange.com/](https://mathematica.stackexchange.com/) and the answer provided by [`user64494`](https://mathematica.stackexchange.com/users/7152/user64494). 

Suppose, we are given two operators $A,B$, and we know their [scale relative graphs](https://arxiv.org/pdf/1902.09788.pdf), denoted by $\mathcal{G}(A)$ and $\mathcal{G}(B)$, respectively. We are interested to figure out the scaled relative graph of their composition $AB:x\mapsto A(B(x))$. For simplicity, we will assume that the regularity conditions under which $\mathcal{G}(AB)=\mathcal{G}(A)\mathcal{G}(B)$ hold; for more details, please see the paper on [scaled relative graph by Ryu et al](https://arxiv.org/pdf/1902.09788.pdf).  

---

Table of contents
\toc

----


### Setup

As an example, let us assume that $A$ is $\beta$-cocoercive and $B$ is $\theta$-averaged ($\theta\in (0,1)$). Then, the scaled relative graph (SRG) of $A$ will be
$$
\mathcal{G}(A)=\left\{z \in \mathbf{C} \mid	\rm{Re}(z) \geq \beta |z|^2\right\},
$$
and the SRG of $B$ will be 
$$
\mathcal G(B)= \left\{ z\in \mathbf{C} \mid 2(1-\theta)\rm{Re}(z) \geq |z|^2 + (1-2\theta)\right\}.
$$
For this example, let us assume, $\beta = 1/2, \theta = 1/2$. 

### Drawing individual SRGs

Let us first draw the individual SRGs for $A$ and $B$. 

**SRG of $A$.** First, we draw the SRG of $A$ by typing the following code. We denote the real and imaginary variables corresponding to $\mathcal G(A)$ by $x,y$, respectively. 

```mathematica
\[Beta] = 1/2;

(*Inequality that defines the SRG of A*)
gAineq = ComplexExpand[Re[z] - \[Beta]*Abs[z]^2 /. z -> x + I*y]

(*Plot the SRG of A*)
srgA = RegionPlot[gAineq >= 0, {x, -2, 2}, {y, -2, 2}]
```

![image-20201118084023741](https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/_assets/image-20201118084023741.png)

We see that 
$$
\mathcal G(A) = \left\{ (x,y) \mid \texttt{gAineq}(x,y) \geq 0\right\}.
$$


**SRG of $B$.** First, we draw the SRG of $B$ by typing the following code. We denote the real and imaginary variables corresponding to $\mathcal G(B)$ by $s,t$, respectively. 

```mathematica 
\[Theta] = 1/2;

(*Inequality that defines the SRG of B*) 
gBineq = ComplexExpand[
  2 (1 - \[Theta]) Re[w] - Abs[w]^2 - (1 - 2 \[Theta]) /. w -> s + I t]

(*Plot the SRG of B*)
srgB = RegionPlot[gBineq >= 0, {s, -2, 2}, {t, -2, 2}]
```

![image-20201118084104555](https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/_assets/image-20201118084104555.png)

Similarly, 
$$
\mathcal G(B) = \left\{ (s,t) \mid \texttt{gBineq}(s,t) \geq 0\right\}.
$$

### Quantifier definition of $\mathcal {G}(A) \mathcal {G}(B)$

We now write down the quantifier definition of $\mathcal {G}(A) \mathcal {G}(B)$. 

By definition: 
$$
z=x+ i y \in \mathcal {G}(A), w=s+it\in \mathcal{G}(B)\Leftrightarrow zw = u + iv \in \mathcal {G}(A) \mathcal {G}(B).
$$


We can break down the equivalence above more by writing it down explicitly into real and imaginary parts.

```mathematica
(*Find out what zw is*)
 zw = ComplexExpand[(x + I y) (s + I t), Reals]

(*Find the imaginary component of zw*) 
v = Plus @@ (Cases[zw, _Complex _]/I)

(*Find the real component of zw*) u = Expand[zw - I v]
```

<img src="https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/_assets/image-20201118084143443.png" alt="image-20201118084143443" style="zoom:67%;" />

So, in terms of quantifier notation, the set description of $\mathcal {G}(A) \mathcal {G}(B)$ will be the following:


$$
\begin{align*}
 & \mathcal{G}(AB)\\
 & =\mathcal{G}(A)\mathcal{G}(B)\\
 & =\left\{ (u,v)\mid u=sx-ty,v=tx+sy,\texttt{gAineq}(x,y)\geq0,\texttt{gBineq}(s,t)\geq0\right\}. \qquad(1)
\end{align*}
$$



### Finding $\mathcal{G}(AB)$ explicitly

In $(1)$, we have $\mathcal{G}(AB)$ in a parametric form where we do not have $u,v$ explicitly, but it is expressed in terms of $x,y$. To figure out the explicit description of $\mathcal{G}(AB)$, we use quantifier elimination technique in `Mathematica`. First step is to observe that: 


$$
\begin{align*}
 & (u,v)\in\mathcal{G}(AB)\\
\Leftrightarrow & \exists_{x,y,s,t\in\mathbf{R}}\left(u=sx-ty,\;v=tx+sy,\;\texttt{gAineq}(x,y)\geq0,\;\texttt{gBineq}(s,t)\geq0\right).
\end{align*}
$$


We next write down this quantifier definition in `Mathematica`. (The rest of the code is self-contained, by only changing the defining inequalities for $A,B$ we can find  $\mathcal{G}(AB)$ for any $A,B$.)

```mathematica
(* Clear the memory *)
ClearAll["Global`*"];

(* Set value of \[Beta] and \[Theta] *)
\[Beta] = 1/2; \[Theta] = 1/2;

(* Define gAineq[x,y], which defines the SRG of operator A via gAineq[x,y] >= 0 *)
gAineq[x_, y_] := 
  ComplexExpand[Re[z] - \[Beta]*Abs[z]^2 /. z -> x + I*y];

(* Define gBineq[s,t], which defines the SRG of operator B via gBineq[s,t] >= 0 *)
gBineq[s_, t_] := 
  ComplexExpand[
   2 (1 - \[Theta]) Re[w] - Abs[w]^2 - (1 - 2 \[Theta]) /. 
    w -> s + I t];

(* Define the quantifier definition for G(AB) *)
quantgAB = 
 Exists[{x, y, s, t}, 
  u == x s - y t && v == x t + y s && gAineq[x, y] >= 0 && 
   gBineq[s, t] >= 0]
```

which gives the output: 

<img src="https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/_assets/image-20201118084216189.png" alt="image-20201118084216189" style="zoom: 50%;" />

Next, we can find the explicit form of $\mathcal{G}(AB)$ in $(u,v)$ by using the `Resolve` command.

```mathematica
(* This will find the explicit form of SRG of AB *)
gAB = Resolve[quantgAB, Reals]
```

which produces the output: 

![image-20201118084318543](https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/_assets/image-20201118084318543.png)

Finally, we can plot $\mathcal{G}(AB)$ as follows.

```mathematica
Region[ImplicitRegion[gAB, {u, v}]]
```

The output is: 

![image-20201118084340687](https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/_assets/image-20201118084340687.png)



