@def title = "An improved lower-complexity bound for the class of convex and Lipschitz continuous functions"
@def published ="November 11, 2021"
@def tags =["optimization"]


# An improved lower-complexity bound for the class of convex and Lipschitz continuous functions

**Shuvomoy Das Gupta**

*November 11, 2021*

In this blog, we discuss how to derive a lower-complexity bound for the class of convex and Lipschitz continuous functions based on [1], which improves the previous lower bound I showed in [this blog](https://shuvomoy.github.io/blogs/posts/A-lower-complexity-bound-for-the-class-of-convex-and-Lipschitz-continuous-functions/). 

**Notation used in the blog.** We use the notation $[1:N]=\{1,2,3,\ldots,N\}$. For a set $S\in \mathbf{R}^d$, its convex hull is denoted by $\mathbf{conv}S$.   For a $d$-dimenstional vector  
$$x = \{x[1],x[2],\ldots,x[d]\} \in \mathbf{R}^d,$$  
define the index set (elements are sorted from smallest to largest):
$$I_{[1:N+1]}(x)=\{\ell\in[1:N+1]\mid x[\ell]=\max\{x[1],x[2],\ldots x[N+1]\}\},$$ *i.e.*, any element of $I_{[1:N+1]}(x)$ corresponds to an index of a maximal component of vector $x$ searched over its first $N+1$ components. The first element in $I_{[1:N+1]}(x)$ is defined by 
$$i_{x}=\min I_{[1:N+1]}(x).$$ 
Finally, For any $i\in [1:d]$, $e_i \in \mathbf{R}^d$ corresponds to the unit vector that has 1 in its $i$-th coordinate and 0 everywhere else.


---

**Table of contents**

\toc

---

#### Large-scale convex optimization problem.

Consider an large-scale convex optimization problem of the form $$\begin{array}{ll}
\textrm{minimize} & f(x),\quad(P)\end{array}$$ where $x\in\mathbf{R}^{d}$ is the decision variable and $d$ is a very large number which is common in today's machine learning problems. The function $f:\mathbf{R}^{d}\to\mathbf{R}$ is assumed to be convex with its subgradient norm bounded by some constant $L>0$, and $x_{\star}$ is an optimal solution to $(P)$.

#### The concept of oracle.

We assume that information on this objective function $f$ is gained through a first-order oracle. Given an input point $x\in\mathbf{R}^{d},$this oracle returns

-   $f(x)$, which is the value of the objective function, and

-   $f^{\prime}(x)$, which is one element in its subdifferential set $\partial f(x).$

#### Algorithm in consideration to solve (P).

Suppose we want to find a solution to $(P),$ so we have to run some iterative algorithm initialized at some point $x_{0}\in\mathbf{R}^{d}$ and, then, using some algorithmic rules, compute $x_{1},x_{2},\ldots,x_{N}$, *i.e.*, $N$ additional iterates. Because we have a large-scale setup, it is safe to assume that $N\leq d-1$ (if the decision variable is a vector with 10 billion components, it is not realistic to compute $10$ billion iterates of any algorithm!). What type of algorithm should we consider?

To solve this optimization problem, we consider algorithm that satisfies the condition $$\left(\forall i\in[1:N]\right)\quad x_{i}\in x_{0}+\textrm{span}\left\{ f^{\prime}(x_{0}),\ldots,f^{\prime}(x_{i-1})\right\} ,\quad\textrm{(Span-Cond)}$$ where $g_{i}\in\partial f(x_{i})$. This condition holds for majority of practical methods.

#### A lower bound result. 

For any $L,R>0$, $d\in\mathbf{N}$ with $N\leq d-1$, and any starting point $x_{0}\in\mathbf{R}^{d}$, there exist (i) a function $f:\mathbf{R}^{d}\to\mathbf{R}$, which is convex and $L$-Lipschitz continuous, and (ii) a first order oracle providing subgradient and function value information such that we have the lower bound on the objective inaccuracy

$$f(x_{N})-f(x_{\star})\geq\frac{LR}{\sqrt{N+1}},$$ for any first-order method satisfying (Span-Cond).

#### Proof of the lower bound result.

The proof proceeds by constructing a "worst-case" function, on which any first-order method satisfying (Span-Cond) will not be able to improve its initial objective value $f(x_{0})$ during the first $N$ computed iterations.

#### Define the worst case function $f$.

First, we define two convex functions $g$​ and $h$​ as follows. Define the convex function $$g(x)\coloneqq g_{N+1}(x)=\max\{x[1],x[2],\ldots x[N+1]\},$$​ where $x[i]$​ denotes the $i$​-th component of the vector $x$​. Using subdifferential calculus [2], we can write down the closed-form expression of $g$​'s subdifferential at a point $x\in\mathbf{R}^{d}$​, given by $$\begin{aligned}
\partial g(x) & =\mathbf{conv}\{e_{k}\mid k\in I_{[1:N+1]}(x)\}\\
 & =\left\{ \sum_{k\in I_{[1:N+1]}(x)}\lambda_{k}e_{k}\mid\left(\forall k\in I_{[1:N+1]}(x)\right)\;\lambda_{k}\geq0,\sum_{k\in I(x)}\lambda_{k}=1\right\} \end{aligned}$$​ where (recall from notation) $$I_{[1:N+1]}(x)=\{\ell\in[1:N+1]\mid x[\ell]=\max\{x[1],x[2],\ldots x[N+1]\}\},$$​ *i.e.*, any element of $I_{[1:N+1]}(x)$​ corresponds to an index of a maximal component of vector $x=\{x[1],\ldots,x[N+1],\ldots,x[d]\}$​ searched over its first $N+1$​ components. Next, define the function $$h(x)\coloneqq h_{N+1}(x)=\|x\|_{2}-R\left(1+\frac{1}{\sqrt{N+1}}\right),$$​ which has subdifferential [2] $$\partial h(x)=\begin{cases}
\frac{x}{\|x\|_{2}}, & \textrm{if }x\neq0,\\
\{y\mid\|y\|_{2}\leq1\}, & \textrm{if }x=0.
\end{cases}$$​ Finally define a function $f$​, which is going to be our worst-case function: $$f(x)\coloneqq f_{N+1}(x)=L\max\left\{ g(x),h(x)\right\} .$$​ Note that $f$​ is a convex function because it is point-wise maximum of two convex functions. Using subdifferential calculus, we can write down the closed-form expression of $f$​'s subdifferential at a point $x\in\mathbf{R}^{d}$​, given by [2] : 

$$
\partial f(x)= L \times \begin{cases}
 \partial g(x), & \textrm{if } g(x) > h(x),\\
 \partial h(x), & \textrm{if } g(x) < h(x),\\
 \mathbf{conv}\{g^{\prime}(x),h^{\prime}(x)\mid g^{\prime}(x)\in\partial g(x),h^{\prime}(x)\in\partial h(x)\}, & \textrm{if } g(x) = h(x).
\end{cases}
$$


#### A handy inequality to bound the subgradient of $f$.

Hence, if we consider a point $x$ and any subgradient of $g$ at $x,$ denoted by $g^{\prime}(x)$, it can be expressed as $$g^{\prime}(x)=\sum_{k\in I_{[1:N+1]}(x)}\lambda_{k}e_{k}:\textrm{ where }\left(\forall k\in I_{[1:N+1]}(x)\right)\;\lambda_{k}\geq0,\sum_{k\in I(x)}\lambda_{k}=1,$$ with the norm satisfying: $$\begin{aligned}
\|g^{\prime}(x)\|_{2} & =\|\sum_{k\in I_{[1:N+1]}(x)}\lambda_{k}e_{k}\|_{2}\\
 & \overset{a)}{\leq}\|\sum_{k\in I_{[1:N+1]}(x)}\lambda_{k}e_{k}\|_{2}\\
 & \overset{b)}{\leq}\sum_{k\in I_{[1:N+1]}(x)}\lambda_{k}\underbrace{\|e_{k}\|_{2}}_{=1}\\
 & =\underbrace{\sum_{k\in I_{[1:N+1]}(x)}\lambda_{k}}_{=1}\\
 & =1\quad(1)\end{aligned}$$ where $a)$ and $b)$ use the triangle inequality and the fact that $\lambda_{k}\geq0$ for all $k\in I_{[1:N+1]}(x)$. Also, it is clear from the subdifferential of $h$ that, any $h^{\prime}(x)\in\partial h(x)$ would satisfy $$\|h^{\prime}(x)\|_{2}\leq1.\quad(2)$$ Hence, from the subdifferential of $f$ above, and (1), (2), we find that any $f^{\prime}(x)\in\partial f(x)$ will satisfy: $$\|f^{\prime}(x)\|_{2}\leq L.\quad(3)$$

#### Verifying that $f$ is $L$-Lipschitz continuous.

Let us now verify that $f$ is indeed Lipschitz continuous. For any $x,y$ and any $f^{\prime}(x)\in\partial f(x)$, from the subgradient inequality, we have $$\begin{aligned}
f(y) & \geq f(x)+\left\langle f^{\prime}(x)\mid y-x\right\rangle \\
\Leftrightarrow f(x)-f(y) & \leq-\left\langle f^{\prime}(x)\mid y-x\right\rangle \\
 & =\left\langle f^{\prime}(x)\mid x-y\right\rangle \\
 & \overset{a)}{\leq}\|f^{\prime}(x)\|_{2}\|x-y\|_{2}\\
 & \leq L\|x-y\|_{2}\end{aligned}$$ where $a)$ uses (3). The last inequality is symmetric with respect to $x$ and $y$, hence for any $x,y$, we also have: $$f(y)-f(x)\leq L\|y-x\|_{2}.$$ Thus, combining the last two inequalities, for any $x,y$, we have $$|f(y)-f(x)|=\max\{f(y)-f(x),-\left(f(y)-f(x)\right)\}\leq L\|y-x\|_{2},$$ *i.e.*, $f$ is $L$-Lipschitz continuous.

#### Defining the oracle.

Now let us define the oracle. When asked for a first-order information about the function at a point $x$​, the oracle returns the function value $f(x)$​, and one subgradient of $f$​ at $x$​ given by (note that it lies in $\partial f(x)$​ defined above) 

$$
f^{\prime}(x)=L \times \begin{cases}
e_{i_{x}}, & \textrm{if } g(x) \geq h(x),\\
x/\|x\|_{2}, & \textrm{if } g(x) < h(x),
\end{cases}
$$

where (recall from notation) $i_{x}$  is the *first* coordinate of $x$ that satisfies $x[i_{x}]=\max_{j\in[1:N+1]}x[j],$ *i.e.*, 

$$i_{x}=\min I_{[1:N+1]}(x).$$

We will see soon that this oracle is providing us the worst possible subgradient for reducing function value in $N$  iterates, and this type of oracle is called *resisting oracle.*

#### Constructing a global minimum of $f$.

Define the point $x_{\star}$ as: $$x_{\star}=\begin{cases}
-\frac{R}{\sqrt{N+1}}, & \textrm{if }i\in[1:N+1]\\
0, & \textrm{else}.
\end{cases}\qquad\textrm{(OptPnt)}$$ which has norm $$\|x_{\star}\|_{2}=R,\qquad\textrm{(NrmOptPnt)}$$ along with $$g(x_{\star})=-\frac{R}{\sqrt{N+1}}$$ and $$h(x_{\star})=\|x_{\star}\|_{2}-R\left(1+\frac{1}{\sqrt{N+1}}\right)=-\frac{R}{\sqrt{N+1}}=g(x_{\star}).$$ As a result, $$f(x_{\star})=L\times\max\{g(x_{\star}),h(x_{\star})\}=-\frac{LR}{\sqrt{N+1}}.$$ We claim that $x_{\star}$ is a global minimum of $f$. First, observe that: $$\begin{aligned}
I_{[1:N+1]}(x_{\star}) & =\{\ell\in[1:N+1]\mid x_{\star}[\ell]=\max\{x_{\star}[1],x_{\star}[2],\ldots x_{\star}[N+1]\}\}\\
 & =\{1,2,\ldots,N+1\},\end{aligned}$$ as all the first $N+1$ components of $N$ are the same. Next note that, at $x_{\star}$, $f$ has the subdifferential $$\begin{aligned}
 & \partial f(x_{\star})\\
 & =L\times\mathbf{conv}\{g^{\prime}(x_{\star}),h^{\prime}(x_{\star})\mid g^{\prime}(x_{\star})\in\partial g(x_{\star}),h^{\prime}(x_{\star})\in\partial h(x_{\star})\}\\
 & =L\times\{\alpha g^{\prime}(x_{\star})+(1-\alpha)h^{\prime}(x_{\star})\mid\alpha\in[0,1],g^{\prime}(x_{\star})\in\partial g(x_{\star}),h^{\prime}(x_{\star})\in\partial h(x_{\star})\}\\
 & =L\times\left\{ \alpha\sum_{k\in I_{[1:N+1]}(x)}\lambda_{k}e_{k}+(1-\alpha)\frac{x_{\star}}{\|x_{\star}\|_{2}}\mid\alpha\in[0,1],\left(\forall k\in I_{[1:N+1]}(x)\right)\;\lambda_{k}\geq0,\sum_{k\in I(x)}\lambda_{k}=1\right\} .\end{aligned}$$ In the last equation, if we set $\alpha\coloneqq\left((1/\sqrt{N+1})+1\right)^{-1}$ and $\lambda_{k}\coloneqq\lambda=1/(N+1)$ for all $k\in[1:N+1],$ then we have one subdifferential of $\partial f(x_{\star})$ as follows: $$\begin{aligned}
\partial f(x_{\star}) & \ni L\times(\alpha\sum_{k\in I_{[1:N+1]}(x)}\lambda_{k}e_{k}+(1-\alpha)\frac{x_{\star}}{\|x_{\star}\|_{2}})\\
 & =L\alpha\lambda\left(\underbrace{1}_{\textrm{index}[1]},\underbrace{1}_{\textrm{index}[2]},\ldots,\underbrace{1}_{\textrm{index}[N+1]},0,\ldots,\underbrace{0}_{\textrm{index}[d]}\right)+\\
 & L(1-\alpha)\left(\underbrace{-\frac{1}{\sqrt{N+1}}}_{\textrm{index}[1]},\underbrace{-\frac{1}{\sqrt{N+1}}}_{\textrm{index}[2]},\ldots,\underbrace{-\frac{1}{\sqrt{N+1}}}_{\textrm{index}[N+1]},0,\ldots,\underbrace{0}_{\textrm{index}[d]}\right)\\
 & =L\alpha\frac{1}{N+1}\left(\underbrace{1}_{\textrm{index}[1]},\underbrace{1}_{\textrm{index}[2]},\ldots,\underbrace{1}_{\textrm{index}[N+1]},0,\ldots,\underbrace{0}_{\textrm{index}[d]}\right)+\\
 & L(1-\alpha)\frac{1}{\sqrt{N+1}}\left(\underbrace{-1}_{\textrm{index}[1]},\underbrace{-1}_{\textrm{index}[2]},\ldots,\underbrace{-1}_{\textrm{index}[N+1]},0,\ldots,\underbrace{0}_{\textrm{index}[d]}\right)\\
 & =(0,0,\ldots,0),\end{aligned}$$ where in the last line we have used the value of $\alpha$.

```mathematica
(*Mathematica code*)
alpha = (1 + 1/Sqrt[n + 1])^-1; 
alpha/(n + 1) - (1 - alpha)/Sqrt[n + 1] // Simplify
```

#### Performance of the oracle in reducing function value. 

**Computing iterate $x_1$.** Let us initialize our algorithm with the initial value $x_{0}=0,$ which satisfies $\|x_{\star}-x_{0}\|_{2}=R$. For this $x_{0},$ we have $$g(x_{0})=\max\{0,\ldots,0\}=0\geq h(x_{0})=\|\mathbf{0}\|_{2}-R\left(1+\frac{1}{\sqrt{N+1}}\right)=-R\left(1+\frac{1}{\sqrt{N+1}}\right),$$ and $$\begin{aligned}
I_{[1:N+1]}(x_{0}) & =\{\ell\in[1:N+1]\mid x[\ell]=\max\{0,0,\ldots,0\}\\
 & =\{1,2,\ldots,N+1\},\\
i_{x_{0}} & =\min I(x_{0})=1.\end{aligned}$$

We ask the resisting oracle about first-order information about the function at $x_{0}$​. It returns $$\begin{aligned}
f(x_{0}) & =L\max\left\{ g(x_{0}),h(x_{0})\right\} \\
 & =L\times0=0\end{aligned}$$​ and $$f^{\prime}(x_{0})=L\times e_{i_{x_{0}}}=Le_{1}.$$​ So, $x_{1}$​, which is the first computed iterate by our first-order method is going to satisfy (Span-Cond): 

$$\begin{aligned}
x_{1} & \in x_{0}+\textrm{span}\left\{ f^{\prime}(x_{0})\right\} \\
 & =0+\textrm{span}\{Le_{1}\}\\
 & =\textrm{span}\{e_{1}\}\\
\Rightarrow x_{1} & =\xi_{1}e_{1}\textrm{ for some }\xi_{1}\in\mathbf{R}.\end{aligned}$$​

For the point $x_{1},$​ we have 

$$g(x_{1})=\max\{\xi_{1}e_{1},\ldots,0\}=\begin{cases}
\xi_{1}, & \textrm{if }\xi_{1}\geq0\\
0, & \textrm{if }\xi_{1}<0
\end{cases}$$​​ 

and 

$$h(x_{1})=\|\xi_{1}e_{1}\|_{2}-R\left(1+\frac{1}{\sqrt{N+1}}\right)=\vert\xi_{1}\vert-R\left(1+\frac{1}{\sqrt{N+1}}\right).$$

**Computing iterate $x_{2}$​.** Consider the case $g(x_{1})\geq h(x_{1}),$​ then 

$$\begin{aligned}
I_{[1:N+1]}(x_{1}) & =\{\ell\in[1:N+1]\mid x[\ell]=\max\{\xi_{1},0,\ldots,0\}\}\\
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
\end{cases}\end{aligned}$$​​

Now we ask the oracle about first-order information about the function at $x_{1}$. It returns

$$\begin{aligned}
f(x_{1}) & =L\times\max\{x_{1}[1],x_{1}[2],\ldots,x_{1}[N+1]\}\\
 & \geq L\times x_{1}[N+1],\end{aligned}$$​ 
 and 
 $$\begin{aligned}
f^{\prime}(x_{1}) & =L\times e_{i_{x_{1}}}\\
 & =L(e_{1}\textrm{ or }e_{2})\end{aligned}$$​ 
 Hence, for the case $g(x_{1})\geq h(x_{1}),$​ the second iterate $x_{2}$​ is going to satisfy (Span-Cond): $$\begin{aligned}
x_{2} & \in x_{0}+\textrm{span}\{f^{\prime}(x_{0}),f^{\prime}(x_{1})\}\\
 & =\textrm{span}\{Le_{1},L(e_{1}\textrm{ or }e_{2})\}\\
 & =\textrm{span}\{e_{1},e_{2}\}.\end{aligned}$$​ 

Now consider the case $g(x_1)<h(x_1)$then 

$$\begin{aligned}
f(x_{1}) & =\vert\xi_{1}\vert-R\left(1+\frac{1}{\sqrt{N+1}}\right)\\
 & =L\times\max\{g(x_{1}),h(x_{1})\}\textrm{ (by definition)}\\
 & \geq L\times g(x_{1})\textrm{ (by definition)}\\
 & =L\times\max\{x_{1}[1],x_{1}[2],\ldots,x_{1}[N+1]\}\\
 & \geq L\times x_{1}[N+1],\end{aligned}$$​​ 

and

 $$f^{\prime}(x_{1})=Lx_{1}/\|x_{1}\|_{2}=L\xi_{1}e_{1}/\vert\xi_{1}\vert=(L\xi_{1}/\vert\xi_{1}\vert)e_{1}.$$​​ 

So for the case $g(x_1) < h(x_1)$​​​,​​​ we have $$x_{2}=\textrm{span}\{e_{1}\}.$$​​​​​​ So, irrespective of $g(x_{1})\geq h(x_{1})$​​​​​​ or $g(x_1)<h(x_1),$​​​​​​ we have $$x_{2}=\textrm{span}\{e_{1},e_{2}\}.$$​​​​​​

**Generalization about the iterates.** Continuing in this manner, we can show that for any $i\in[1:N],$ we have $$x_{i}\in\textrm{span}\{e_{1},e_{2},\ldots,e_{i}\},$$ and as a result, components of $x_{i}$ corresponding to indices $N+1,\ldots,d$ will be zero for all $i\in[1:N]$, *i.e.*, $$x_{i}[N+1]=\ldots=x_{i}[d]=0.$$ Hence, the objective value of the iterates $x_{1},\ldots,x_{N}$ is going to satisfy for all $i\in[1:N]$ $$\begin{aligned}
f(x_{i}) & =L\times\max\{g(x_{i}),h(x_{i})\}\textrm{ (by definition)}\\
 & \geq L\times g(x_{i})\textrm{ (by definition)}\\
 & =L\times\max\{x_{i}[1],x_{i}[2],\ldots,x_{i}[N+1]\}\\
 & \geq L\times x_{i}[N+1]=0,\qquad\textrm{(ObjNotImprv)}\end{aligned}$$ where in the last line we have used the simple fact that the maximum element of a list will be greater than any element of the list. This is interesting: because it says that no matter what we do, we cannot improve the initial objective value $f(x_{0})=0$ for the iterates $x_{1},\ldots,x_{N}$! 

 **Final lower bound.** Finally, applying (ObjNotImprv) for the $N$-th iterate, we get $$\begin{aligned}
f(x_{N})-f(x_{\star}) & \geq0-f(x_{\star})\\
 & =-\left(-\frac{LR}{\sqrt{N+1}}\right)\\
 & =\frac{LR}{\sqrt{N+1}}\qquad\textrm{(ObjGap)}\end{aligned}$$

thus completing our proof.

#### References

[1] Drori, Yoel, and Teboulle, Marc. An optimal variant of Kelley's cutting-plane method. Mathematical Programming 160.1 (2016): 321-351.

[2] Beck, Amir. First-order methods in optimization. Society for Industrial and Applied Mathematics, 2017.
