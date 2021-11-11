@def title = "A lower-complexity bound for the class of convex and Lipschitz continuous functions"
@def published ="November 8, 2021"
@def tags =["optimization"]



# A lower-complexity bound for the class of convex and Lipschitz continuous functions

**Shuvomoy Das Gupta**

*November 8, 2021*

In this blog, we discuss how to derive a lower-complexity bound for the class of convex and Lipschitz continuous functions. This blog is based on [1, Section 3.2.1].
~~~
<label for="mn-demo-s1" class="margin-toggle">&#8853;</label>
<input type="checkbox" id="mn-demo" class="margin-toggle"/>
<span class="marginnote">
[1] Y. Nesterov. Introductory lectures on convex optimization: A basic course.
Kluwer Academic Publishers, 2004.
</span>
~~~

---

**Table of contents**

\toc

---
### Background

#### Large-scale convex optimization problem

Consider a large-scale convex optimization problem of the form 

$$\begin{array}{ll}
\textrm{minimize} & f(x),\quad(P)\end{array}$$ 

where $x\in\mathbf{R}^{d}$ is the decision variable and $d$ is a very large number which is common in today's machine learning problems. The function $f$ is assumed to be convex with its subgradient norm bounded by some constant $L>0$ over a large enough ball of radius $R$ around $x_{\star}$, where $x_{\star}$ is an optimal solution to $(P)$.

#### The concept of oracle

We assume that information on this objective function $f$ is gained through a first-order oracle. Given an input point $x\in\mathbf{R}^{d},$ this oracle returns

-   $f(x)$, which is the value of the objective function, and

-   $f^{\prime}(x)$, which is one element in its subdifferential set $\partial f(x).$

#### Algorithm in consideration to solve (P)

Suppose we want to find a solution to $(P),$ so we have to run some iterative algorithm initialized at some point $x_{0}\in\mathbf{R}^{d}$ and, then, using some algorithmic rules, compute $x_{1},x_{2},\ldots,x_{N}$, *i.e.*, $N$ additional iterates. Because we have a large-scale setup, it is safe to assume that $N\leq d-1$ (if the decision variable is a vector with 10 billion components, it is not realistic to compute $10$ billion iterates of any algorithm!). What type of algorithm should we consider?

To solve this optimization problem, we consider algorithm that satisfies the condition 
~~~
<label for="mn-demo-s2" class="margin-toggle">&#8853;</label>
<input type="checkbox" id="mn-demo" class="margin-toggle"/>
<span class="marginnote">
We use the notation \([1:N]={1,2,\ldots,N}\).
</span>
~~~

 $$\left(\forall i\in[1:N]\right)\quad x_{i}\in x_{0}+\textrm{span}\left\{ f^{\prime}(x_{0}),\ldots,f^{\prime}(x_{i-1})\right\} ,\quad\textrm{(Span-Cond)}$$ 

where $g_{i}\in\partial f(x_{i})$. This condition holds for majority of practical methods.

We next present the main lower bound result, and then we prove it.

---
### A lower bound result

For any $L,R>0$​, $d\in\mathbf{N}$​ with $N\leq d-1$​, and any starting point $x_{0}\in\mathbf{R}^{d}$​, there exist (i) a function $f:\mathbf{R}^{d}\to\mathbf{R}$​, which is convex and $L$​-Lipschitz continuous on the ball $B(x_{\star};R)=\{x\in\mathbf{R}^{d}\mid\|x-x_{\star}\|_{2}^{2}\leq R^{2}\}$​ and (ii) a first order oracle providing subgradient and function value information such that we have the lower bound on the objective inaccuracy

 $$f(x_{N})-f(x_{\star})\geq\frac{LR}{2(2+\sqrt{N+1})},$$ 

 for any first-order method satisfying (Span-Cond).

---


### Proof of the lower bound result

The proof proceeds by constructing a "worst-case" function, on which any first-order method satisfying (Span-Cond) will not be able to improve its initial objective value $f(x_{0})$ during the first $N$ computed iterations.

#### Define the worst case function $f$

Define the function: 

$$f(x) := \frac{\mu}{2}\|x\|_{2}^{2}+\gamma\max\{x[1],x[2],\ldots, x[N+1]\},$$ 

where $x[i]$ denotes the $i$-th component of the vector $x$, and $\gamma$ and $\mu$ are some positive numbers. For now, we do not specify what $\mu$ and $\gamma$ are, we will automatically find their values in terms of $L$ and $R$ down the line. Note that $f$ is convex by definition. Using subdifferential calculus, we can write down the closed-form expression of $f$'s subdifferential at a point $x\in\mathbf{R}^{d}$, given by 
~~~
<label for="mn-demo-s3" class="margin-toggle">&#8853;</label>
<input type="checkbox" id="mn-demo" class="margin-toggle"/>
<span class="marginnote">
For any \(i\in\{1,2,\ldots\}\), \(e_i\in\mathbf{R}^d\) corresponds to the unit vector that has 1 in its \(i\)-th coordinate and 0 everywhere else.
</span>
~~~

$$\begin{aligned}
\partial f(x) & =\mu x+\gamma\mathbf{conv}\{e_{k}\mid k\in I(x)\}\\
 & =\mu x+\gamma\left\{ \sum_{k\in I(x)}\lambda_{k}e_{k}\mid\left(\forall k\in I(x)\right)\;\lambda_{k}\geq0,\sum_{k\in I(x)}\lambda_{k}=1\right\} \end{aligned}$$

 where $$I(x)=\{\ell\in[1:N+1]\mid x[\ell]=\max\{x[1],x[2],\ldots ,x[N+1]\}\},$$ *i.e.*, any element of $I(x)$ corresponds to an index of a maximal component of vector $x=\{x[1],\ldots,x[N+1],\ldots,x[d]\}$ searched over its first $N+1$ components. 

#### A handy inequality to bound the subgradient of $f$

Hence, if we consider a point $x\in B(x_{\star};R)$ and any subgradient of $f$ at $x,$ denoted by $f^{\prime}(x)$, it can be expressed as $$f^{\prime}(x)=\mu x+\gamma\sum_{k\in I(x)}\lambda_{k}e_{k}:\textrm{ where }\left(\forall k\in I(x)\right)\;\lambda_{k}\geq0,\sum_{k\in I(x)}\lambda_{k}=1,$$ with the norm satisfying: 

$$\begin{aligned}
\|f^{\prime}(x)\|_{2} & =\|\mu x+\gamma\sum_{k\in I(x)}\lambda_{k}e_{k}\|_{2}\\
 & \overset{a)}{\leq}\mu\|x\|_{2}+\gamma\|\sum_{k\in I(x)}\lambda_{k}e_{k}\|_{2}\\
 & \overset{b)}{\leq}\mu\|x\|_{2}+\gamma\sum_{k\in I(x)}\lambda_{k}\underbrace{\|e_{k}\|_{2}}_{=1}\\
 & =\mu\|x\|_{2}+\gamma\underbrace{\sum_{k\in I(x)}\lambda_{k}}_{=1}\\
 & =\mu\|x\|_{2}+\gamma\\
 & \overset{c)}{\leq}\mu\left(R+\|x_{\star}\|_{2}\right)+\gamma\qquad(1)\end{aligned}$$ 

where $a)$ and $b)$ use the triangle inequality and the fact that $\lambda_{k}\geq0$ for all $k\in I(x)$, and $c)$ uses the reverse triangle inequality: 

$$\begin{aligned}
 & \|x\|_{2}-\|x_{\star}\|_{2}\leq\|x-x_{\star}\|_{2}\leq R\\
\Rightarrow & \|x\|_{2}\leq R+\|x_{\star}\|_{2}.\end{aligned}$$

#### Verifying that $f$ is Lipschitz continuous on $B(x_{\star};R)$

Let us now verify that $f$ is indeed Lipschitz continuous on the ball $B(x_{\star};R)=\{x\in\mathbf{R}^{d}\mid\|x-x_{\star}\|_{2}^{2}\leq R^{2}\}$. For any $x,y\in B(x_{\star};R)$ and any $f^{\prime}(x)\in\partial f(x)$, from the subgradient inequality we have 

$$\begin{aligned}
f(y) & \geq f(x)+\left\langle f^{\prime}(x)\mid y-x\right\rangle \\
\Leftrightarrow f(x)-f(y) & \leq-\left\langle f^{\prime}(x)\mid y-x\right\rangle \\
 & =\left\langle f^{\prime}(x)\mid x-y\right\rangle \\
 & \overset{a)}{\leq}\|f^{\prime}(x)\|_2 \|x-y\|_2\\
 & \overset{b)}{\leq}\left(\mu\left(R+\|x_{\star}\|_{2}\right)+\gamma\right)\|x-y\|_2,\end{aligned}$$ 

where $a)$ uses Cauchy--Schwarz inequality and $b)$ uses (1). The last inequality is symmetric with respect to $x$ and $y$, hence for any $x,y\in B(x_{\star};R)$, we also have: 

$$f(y)-f(x)\leq\left(\mu\left(R+\|x_{\star}\|_{2}\right)+\gamma\right)\|y-x\|_2.$$ 

Thus, combining the last two inequalities, for any $x,y\in B(x_{\star};R)$, we have 

$$|f(y)-f(x)|=\max\{f(y)-f(x),-\left(f(y)-f(x)\right)\}\leq\left(\mu\left(R+\|x_{\star}\|_{2}\right)+\gamma\right)\|y-x\|_2,$$ 

*i.e.*, $f$ is $\left(\mu\left(R+\|x_{\star}\|_{2}\right)+\gamma\right)$-Lipschitz continuous on $B(x_{\star};R)$. 

In other words, when we pick the value of $\mu$ and $\gamma$, we need to satisfy the equation: 

$$L=\mu\left(R+\|x_{\star}\|_{2}\right)+\gamma,\qquad\textrm{(ConstraitntOnL)}$$

which we will do in the end.

#### Defining the oracle

Now let us define the oracle. When asked for first-order information about the function at a point $x$, the oracle returns the function value $f(x)$, and one subgradient of $f$ at $x$ given by 

$$f^{\prime}(x)=\mu x+\gamma e_{i_{x}}\in\partial f(x),$$

 where $i_{x}$ is the *first* coordinate of $x$ that satisfies $x[i_{x}]=\max_{j\in[1:N+1]}x[j],$  *i.e.*, 

$$i_{x}=\min I(x).$$ 

We will see soon that this oracle is providing us the worst possible subgradient for reducing function value in $N$ iterates, and this type of oracle is called *resisting oracle.*

#### Constructing a global minimum of $f$

Define the point $x_{\star}$ as: 

$$x_{\star}=\begin{cases}
-\frac{\gamma}{\mu(N+1)}, & \textrm{if }i\in[1:N+1]\\
0, & \textrm{else},
\end{cases}\qquad\textrm{(OptPnt)}$$ 

which has norm 

$$\|x_{\star}\|_{2}=\frac{\gamma}{\mu\sqrt{N+1}}.\qquad\textrm{(NrmOptPnt)}$$ 

Keep in mind that, in the end we need to pick the value of $\gamma$ and $\mu$ in such a way so that the condition 

$$\|x_{0}-x_{\star}\|_{2}\leq R\qquad\textrm{(OptPntDist)}$$ 

is satisfied. We claim that $x_{\star}$ is a global minimum of $f$. First, observe that: 

$$\begin{aligned}
I(x_{\star}) & =\{\ell\in[1:N+1]\mid x_{\star}[\ell]=\max\{x_{\star}[1],x_{\star}[2],\ldots, x_{\star}[N+1]\}\}\\
 & =\{1,2,\ldots,N+1\},\end{aligned}$$ 

as all the first $N+1$ components of $N$ are the same. Next note that, at $x_{\star}$, $f$ has the subdifferential 

$$\partial f(x_{\star})=\mu x_{\star}+\gamma\left\{ \sum_{k\in[1:N+1]}\lambda_{k}e_{k}\mid\left(\forall k\in[1:N+1]\right)\;\lambda_{k}\geq0,\sum_{k\in I(x)}\lambda_{k}=1\right\} .$$

 If we set, $\lambda_{k}=1/(N+1)$ for $i\in[1:N+1],$then one particular element of $\partial f(x_{\star})$ would be: 

$$\begin{aligned}
\partial f(x_{\star}) & \ni\mu x_{\star}+\gamma\overbrace{\sum_{i\in[1:N+1]}\frac{1}{N+1}e_{i}}^{z(\textrm{say})}\\
 & =\mu\left(\underbrace{-\frac{\gamma}{\mu(N+1)}}_{x_{\star}[1]},\underbrace{-\frac{\gamma}{\mu(N+1)}}_{x_{\star}[2]},\ldots,\underbrace{-\frac{\gamma}{\mu(N+1)}}_{x_{\star}[N+1]},0,\ldots,\underbrace{0}_{x_{\star}[d]}\right)+\\
 & \gamma\left(\underbrace{\frac{1}{(N+1)}}_{z[1]},\underbrace{\frac{1}{(N+1)}}_{z[2]},\ldots,\underbrace{\frac{1}{(N+1)}}_{z[N+1]},0,\ldots,\underbrace{0}_{z[d]}\right)\\
 & =(0,0,\ldots,0),\end{aligned}$$ 

hence via Fermat's rule, $x_{\star}$ is a global minimum of $f$. The optimal value of $f$ at $x_{\star}$ is: 

$$\begin{aligned}
f(x_{\star}) & =\frac{\mu}{2}\|x_{\star}\|_{2}^{2}+\gamma\max\{x_{\star}[1],x_{\star}[2],\ldots x_{\star}[N+1]\}\\
 & =\frac{\mu}{2}\frac{\gamma^{2}}{\mu^{2}(N+1)}+\gamma\left(-\frac{\gamma}{\mu(N+1)}\right)\\
 & =\frac{\gamma^{2}}{2\mu(N+1)}-\frac{\gamma^{2}}{\mu(N+1)}\\
 & =-\frac{\gamma^{2}}{2\mu(N+1)}.\end{aligned}$$ 

#### Performance of the oracle in reducing function value

Let us initialize our algorithm with the initial value $x_{0}=0,$ for which we have

 $$\begin{aligned}
I(x_{0}) & =\{\ell\in[1:N+1]\mid x[\ell]=\max\{0,0,\ldots,0\}\\
 & =\{1,2,\ldots,N+1\},\\
i_{x_{0}} & =\min I(x_{0})=1.\end{aligned}$$

We ask the resisting oracle about first-order information about the function at $x_{0}$. It returns 

$$\begin{aligned}
f(x_{0}) & =\frac{\mu}{2}\|x_{0}\|_{2}^{2}+\gamma\max\{x_{0}[1],x_{0}[2],\ldots, x_{0}[N+1]\}\\
 & =0,\end{aligned}$$ 

and 

$$\begin{aligned}
f^{\prime}(x_{0}) & =\mu x_{0}+\gamma e_{i_{x_{0}}}\\
 & =\gamma e_{1}.\end{aligned}$$ 

So, $x_{1}$, which is the first computed iterate by our first-order method is going to satisfy (Span-Cond): 

$$\begin{aligned}
x_{1} & \in x_{0}+\textrm{span}\left\{ f^{\prime}(x_{0})\right\} \\
 & =0+\textrm{span}\{\gamma e_{1}\}\\
 & =\textrm{span}\{e_{1}\}\\
\Rightarrow x_{1} & =\xi_{1}e_{1}\textrm{ for some }\xi_{1}\in\mathbf{R}.\end{aligned}$$

For the point $x_{1},$ we have 

$$\begin{aligned}
I(x_{1}) & =\{\ell\in[1:N+1]\mid x[\ell]=\max\{\xi_{1},0,\ldots,0\}\}\\
 & =\begin{cases}
\{1,2,\ldots,N+1\}, & \textrm{if }\xi_{1}=0,\\
\{2,3,\ldots,N+1\}, & \textrm{if }\xi_{1}<0,\\
\{1\}, & \textrm{if }\xi_{1}>0,
\end{cases}\end{aligned}$$ 

and 

$$\begin{aligned}
i_{x_{1}} & =\min I(x_{1})\\
 & =\begin{cases}
1, & \textrm{if }\xi_{1}=0,\\
2, & \textrm{if }\xi_{1}<0,\\
1, & \textrm{if }\xi_{1}>0.
\end{cases}\end{aligned}$$

Now we ask the oracle about first-order information about the function at $x_{1}$. It returns

$$\begin{aligned}
f(x_{1}) & =\frac{\mu}{2}\|x_{1}\|_{2}^{2}+\gamma\max\{x_{1}[1],x_{1}[2],\ldots,x_{1}[N+1]\}\end{aligned}$$

 and 

$$\begin{aligned}
f^{\prime}(x_{1}) & =\mu x_{1}+\gamma e_{i_{x_{1}}}\\
 & =\mu\xi_{1}\gamma e_{1}+\gamma(e_{1}\textrm{ or }e_{2})\end{aligned}$$

 Hence, the second iterate $x_{2}$ is going to satisfy (Span-Cond): 

$$\begin{aligned}
x_{2} & \in x_{0}+\textrm{span}\{f^{\prime}(x_{0}),f^{\prime}(x_{1})\}\\
 & =\textrm{span}\{\gamma e_{1},\mu\xi_{1}\gamma e_{1}+\gamma(e_{1}\textrm{ or }e_{2})\}\\
 & =\textrm{span}\{e_{1},e_{2}\}.\end{aligned}$$

Continuing in this manner, we can show that for any $i\in[1:N],$ we have $$x_{i}\in\textrm{span}\{e_{1},e_{2},\ldots,e_{i}\},$$ and as a result, components of $x_{i}$ corresponding to indices $N+1,\ldots,d$ will be zero for all $i\in[1:N]$, *i.e.*, 

$$x_{i}[N+1]=\ldots=x_{i}[d]=0.$$

 Hence, the objective value of the iterates $x_{1},\ldots,x_{N}$ is going to satisfy for all $i\in[1:N]$ 

$$\begin{aligned}
f(x_{i}) & =\frac{\mu}{2}\|x_{i}\|_{2}^{2}+\gamma\max\{x_{i}[1],x_{i}[2],\ldots,x_{i}[N+1]\}\\
 & \geq\gamma\max\{x_{i}[1],x_{i}[2],\ldots,x_{i}[N+1]\}\\
 & \geq \gamma x_{i}[N+1]\\
 & =0,\qquad\textrm{(ObjNotImprv)}\end{aligned}$$​

where in the second last line we have used the simple fact that the maximum element of a list will be greater than any element of the list. This is interesting: because it says that no matter what we do, we cannot improve the initial objective value $f(x_{0})=0$ for the iterates $x_{1},\ldots,x_{N}$! Finally, applying (ObjNotImprv) for the $N$-th iterate, we get 

$$\begin{aligned}
f(x_{N})-f(x_{\star}) & \geq0-f(x_{\star})\\
 & =-\left(-\frac{\gamma^{2}}{2\mu(N+1)}\right)\\
 & =\frac{\gamma^{2}}{2\mu(N+1)}.\qquad\textrm{(ObjGap)}\end{aligned}$$

#### Picking values of $\gamma$ and $\mu$

Now we are left to picking the values of $\gamma$ and $\mu$ in terms of $L$ and $R$. For that we need to be mindful of two constraints that $\gamma$ and $\mu$ need to satisfy:

- From (ConstraitntOnL) :

  $$\begin{aligned}
  L & =\mu\left(R+\|x_{\star}\|_{2}\right)+\gamma\\
   & =\mu\left(R+\frac{\gamma}{\mu\sqrt{N+1}}\right)+\gamma,\end{aligned}$$ and

- From (OptPntDist) :

  $$\|x_{0}-x_{\star}\|_{2}^{2}=\|0-x_{\star}\|_{2}^{2}=\frac{\gamma^{2}}{\mu^{2}(N+1)}\leq R^{2},$$ 

  where we have used (NrmOptPnt).

To make our life easier, we consider equality in the last equation, and solver for $\gamma$ and $\mu$. We can do it by hand, but I am lazy, so I have just used Mathematica.

```julia
(*Mathematica code*)
In[1] := Solve[Reduce[{L == \[Mu] (R + \[Gamma]/(\[Mu] Sqrt[\[CapitalNu] + 
            1])) + \[Gamma], \[Gamma]/(\[Mu] Sqrt[\[CapitalNu] + 1]) == 
     R , R > 0, \[CapitalNu] > -1, 
    L > 0, \[Gamma] > 0, \[Mu] > 0}, {\[Gamma], \[Mu]}, 
   Reals], {\[Gamma], \[Mu]}, Reals] // Simplify
```

The solution is: 

$$\begin{aligned}
\gamma & =\frac{L\sqrt{N+1}}{2+\sqrt{N+1}},\\
\mu & =\frac{L}{2R+R\sqrt{N+1}}.\end{aligned}$$ 

Putting the values of $\gamma$ and $\mu$ in (ObjGap), 
```julia
(*Mathematica code*)
In[2] := \[Gamma] = (L Sqrt[1 + \[CapitalNu]])/(2 + Sqrt[1 + \[CapitalNu]]);
\[Mu] = L/(2 R + R Sqrt[1 + \[CapitalNu]]);
ObjGap = \[Gamma]^2/(2 \[Mu] (\[CapitalNu] + 1)) // FullSimplify
```
we find: 

$$\begin{aligned}
f(x_{N})-f(x_{\star}) & \geq\frac{\gamma^{2}}{2\mu(N+1)}\\
 & =\frac{LR}{4+2\sqrt{N+1}},\end{aligned}$$ 

thus completing our proof.

