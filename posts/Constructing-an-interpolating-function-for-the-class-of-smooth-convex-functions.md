\enabletheorems

@def title = "Constructing an interpolating function for the class of smooth convex functions"
@def published ="November 22, 2021"
@def tags =["optimization", "lower-bounds-in-optimization"]

# Constructing an interpolating function for the class of smooth convex functions

**Shuvomoy Das Gupta**

*November 22, 2021*

In this blog, we study constructing interpolating function for $L$​​-smooth convex function due to Yoel Drori from `[3]`.

**Notation and notions.** All norms are 2-norm in this blog. A function is $L$-smooth convex if and only 
$$\left(\forall x,y\in\mathbf{R}^{d}\right)\quad f(y)\geq f(x)+\langle\nabla f(x)\mid y-x\rangle+\frac{1}{2L}\|\nabla f(x)-\nabla f(y)\|^{2}.$$

On the other hand, a function is $\mu$​-strongly convex if and only if 

$$\left(\forall x,y\in\mathbf{R}^{d}\right)\quad f(y)\geq f(x)+\langle\nabla f(x)\mid y-x\rangle+\frac{\mu}{2}\|f^{\prime}(x)-f^{\prime}(y)\|^{2},$$​ 
where $f^{\prime}(\cdot)$​​ denotes a subgradient of $f$​​ at $(\cdot)$​​. The set of all $L$​​-smooth convex functions on $\mathbf{R}^{d}$​​ is denoted by $\mathcal{F}_{0,L}(\mathbf{R}^{d})$​​, and the set of all $\mu$​​-strongly convex functions is denoted by $\mathcal{F}_{\mu,\infty}(\mathbf{R}^{d})$​​​.


---

**Table of contents**

\toc

---



### Sion's minimax theorem

We are going to use the famous result due to Maurice Sion.

\begin{theorem}{Sion's minimax theorem}{sion-mini-max}
Let $X$ be a compact convex set (over which minimization will be performed) and $Y$ be a convex set in $\mathbf{R}^{d}$ (over which supremum will be computed). Suppose $g:X\times Y\to\mathbf{R}$ satisfies the following properties: (i) $g(x,\cdot)$ is upper-semicontinuous and quasi-concave on $Y$ for all $x\in X$, and (ii) $g(\cdot,y)$ is lower-semicontinuous and quasi-convex on $X$ for all $y\in Y$. Then we have 
$$\sup_{y\in Y}\min_{x\in X}g(x,y)=\min_{x\in X}\sup_{y\in Y}g(x,y).$$
\end{theorem}

### Interpolable function

First, we start with the definition of an interpolable function.

\begin{definition}{Interpolable function.}{interpol-fun}
Suppose we are given the set of triplets $\{x_{i},g_{i},f_{i}\}_{i\in I}$ where $I$ is a finite index set and $x_{i}\in\mathbf{R}^{d},g_{i}\in\mathbf{R}^{d},$ and $f_{i}\in\mathbf{R}.$ Then the set $\{x_{i},g_{i},f_{i}\}_{i\in I}$ is $\mathcal{F}_{0,L}(\mathbf{R}^{d})$-interpolable if and only if there exists a function $f\in\mathcal{F}_{0,L}(\mathbf{R}^{d})$ such that for all $i\in I$ we have $f_{i}=f(x_{i})$ and $g_{i}=\nabla f(x_{i})$​.
\end{definition}

### Main result (due to Yoel Drori)

Next, we prove the following result due to Yoel Drori from `[3]`.

\begin{theorem}{Interpolation of smooth convex functions.}{intpl-smth-convex-func}
If $\{x_{i},g_{i},f_{i}\}_{i\in I}\subseteq\mathbf{R}^{d}\times\mathbf{R}^{d}\times\mathbf{R}$​ is $\mathcal{F}_{0,L}(\mathbf{R}^{d})$​-interpolable, then one interpolation function $f\in\mathcal{F}_{0,L}(\mathbf{R}^{d})$​ that can be constructed from $\{x_{i},g_{i},f_{i}\}_{i\in I}$​ is: $$f(x)=\min_{\alpha\in\Delta}\left[\frac{L}{2}\|x-\sum_{i\in I}\alpha_{i}\left(x_{i}-\frac{1}{L}g_{i}\right)\|^{2}+\sum_{i\in I}\alpha_{i}\left(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\right)\right],$$​ where $\Delta=\{\bar{\alpha}\in\mathbf{R}^{\vert I\vert}\mid\bar{\alpha}\geq0,\sum_{i\in I}\bar{\alpha}_{i}=1\}.$​
\end{theorem}

### Proof

The set $\{x_{i},g_{i},f_{i}\}_{i\in I}$ is $\mathcal{F}_{0,L}(\mathbf{R}^{d})$-interpolable if and only if $\{\bar{x}_{i},\bar{g}_{i},\bar{f}_{i}\}_{i\in I}\coloneqq\{g_{i},x_{i},\left\langle x_{i}\mid g_{i}\right\rangle -f_{i}\}_{i\in I}$ is $\mathcal{F}_{1/L,\infty}(\mathbf{R}^{d})$, which is from `[1, Lemma 3.7]`. Also, if $\{\bar{x}_{i},\bar{g}_{i},\bar{f}_{i}\}_{i\in I}$ is $\mathcal{F}_{\mu,\infty}(\mathbf{R}^{d})$-interpolable then one such interpolation function $h\in\mathcal{F}_{\mu,\infty}(\mathbf{R}^{d})$ would be: $$\begin{aligned}
h(\bar{x}) & =\max_{i\in I}\{\bar{f}_{i}+\left\langle \bar{g}_{i}\mid\bar{x}-\bar{x}_{i}\right\rangle +\frac{\mu}{2}\|\bar{x}-\bar{x}_{i}\|^{2}\}\end{aligned}$$ where we have used the fact that $\max_{i\in I}a_{i}=\max_{\alpha\in\Delta}\sum_{i\in I}\alpha_{i}a_{i}.$

So, if $\{\bar{x}_{i},\bar{g}_{i},\bar{f}_{i}\}_{i\in I}\coloneqq\{g_{i},x_{i},\left\langle x_{i}\mid g_{i}\right\rangle -f_{i}\}_{i\in I}$ is $\mathcal{F}_{1/L,\infty}(\mathbf{R}^{d})$-interpolable, then one such interpolation function $h\in\mathcal{F}_{1/L,\infty}(\mathbf{R}^{d})$ would be $$\begin{aligned}
h(\bar{x}) & =\max_{i\in I}\left\{ \bar{f}_{i}+\left\langle \bar{g}_{i}\mid\bar{x}-\bar{x}_{i}\right\rangle +\frac{1}{2L}\|\bar{x}-\bar{x}_{i}\|^{2}\right\} \\
 & =\max_{i\in I}\left\{ \left\langle x_{i}\mid g_{i}\right\rangle -f_{i}+\left\langle x_{i}\mid\bar{x}-g_{i}\right\rangle +\frac{1}{2L}\|\bar{x}-g_{i}\|^{2}\right\} \\
 & =\max_{i\in I}\left\{ \cancel{\left\langle x_{i}\mid g_{i}\right\rangle }-f_{i}+\left\langle x_{i}\mid\bar{x}\right\rangle -\cancel{\left\langle x_{i}\mid g_{i}\right\rangle }+\frac{1}{2L}\|\bar{x}-g_{i}\|^{2}\right\} \\
 & =\max_{i\in I}\left\{ -f_{i}+\left\langle x_{i}\mid\bar{x}\right\rangle +\frac{1}{2L}\|\bar{x}-g_{i}\|^{2}\right\} \\
 & =-\min_{i\in I}\left\{ f_{i}-\left\langle x_{i}\mid\bar{x}\right\rangle -\frac{1}{2L}\|\bar{x}-g_{i}\|^{2}\right\} ,\quad\quad\quad(1)\end{aligned}$$ where in the last line we have used $\max(\cdot)=-\min(-\cdot).$

But, we are looking for an interpolation function that is in $\mathcal{F}_{0,L}(\mathbf{R}^{d}),$ not $\mathcal{F}_{1/L,\infty}(\mathbf{R}^{d})$. How to go from $\mathcal{F}_{1/L,\infty}(\mathbf{R}^{d})$ to $\mathcal{F}_{0,L}(\mathbf{R}^{d})$? To that goal, we use the following results:

\(i\) $f\in\mathcal{F}_{0,L}(\mathbf{R}^{d})$​ if and only if $f$​'s conjugate function $f^{*}(\cdot)=\sup_{y\in\mathbf{R}^{d}}\left[-f(y)+\left\langle \cdot\mid y\right\rangle \right]\in\mathcal{F}_{1/L,\infty}(\mathbf{R}^{d})$​ `[1, Theorem 2.34]`.

\(ii\) for any lower-semicontinuous, proper, and convex function $f$, we have $f=f^{**}$ `[2, Theorem 12.2]`.

Due to (i) and (ii), to find the interpolation function in $\mathcal{F}_{0,L}(\mathbf{R}^{d}),$​​ all we have to do is to compute the conjugate function of $h$​​. In other words, the desired interpolation function in $\mathcal{F}_{0,L}(\mathbf{R}^{d})$​​ would be 

$$\begin{aligned}
f(x) & =h^{*}(x)\\
 & =\sup_{y\in\mathbf{R}^{d}}\left[-h(y)+\left\langle x\mid y\right\rangle \right]\\
 & \overset{(1)}{=}\sup_{y\in\mathbf{R}^{d}}\left[-\left[-\min_{i\in I}\left\{ f_{i}-\left\langle x_{i}\mid y\right\rangle -\frac{1}{2L}\|y-g_{i}\|^{2}\right\} \right]+\left\langle x\mid y\right\rangle \right]\\
 & =\sup_{y\in\mathbf{R}^{d}}\left[\min_{i\in I}\left\{ f_{i}-\left\langle x_{i}\mid y\right\rangle -\frac{1}{2L}\|y-g_{i}\|^{2}+\left\langle x\mid y\right\rangle \right\} \right]\\
 & =\sup_{y\in\mathbf{R}^{d}}\left[\min_{i\in I}\left\{ f_{i}-\frac{1}{2L}\|y-g_{i}\|^{2}+\left\langle x-x_{i}\mid y\right\rangle \right\} \right]\\
 & =\sup_{y\in\mathbf{R}^{d}}\left[\min_{\alpha\in\Delta}\sum_{i\in I}\alpha_{i}\left\{ f_{i}-\frac{1}{2L}\|y-g_{i}\|^{2}+\left\langle x-x_{i}\mid y\right\rangle \right\} \right],\quad\quad\quad(2)\end{aligned}$$​​ 
where in the last line we have used the fact that $\min_{i\in I}a_{i}=\min_{\alpha\in\Delta}\sum_{i\in I}\alpha_{i}a_{i}$​​​. Now, in (2), denote the inner function by $$p(y,\alpha)=\sum_{i\in I}\alpha_{i}\left\{ f_{i}-\frac{1}{2L}\|y-g_{i}\|^{2}+\left\langle x-x_{i}\mid y\right\rangle \right\} ,$$​​​ which is continuous and concave with respect to the maximizing variable $y$​​​ and continuous and convex (in fact linear) with respect to the minimizing variable $\alpha$​​​. Finally, $\Delta$​​​ the minimizing set is compact and convex, and the $\mathbf{R}^{d}$​​​ the maximizing set is convex. Hence, applying Sion's minimax theorem we have: 

$$\begin{aligned}
 & f(x)\\
= & \sup_{y\in\mathbf{R}^{d}}\left[\min_{\alpha\in\Delta}\sum_{i\in I}\alpha_{i}\left\{ f_{i}-\frac{1}{2L}\|y-g_{i}\|^{2}+\left\langle x-x_{i}\mid y\right\rangle \right\} \right]\\
= & \min_{\alpha\in\Delta}\left[\sup_{y\in\mathbf{R}^{d}}\sum_{i\in I}\alpha_{i}\left\{ f_{i}-\frac{1}{2L}\|y-g_{i}\|^{2}+\left\langle x-x_{i}\mid y\right\rangle \right\} \right].\quad\quad\quad(3)\end{aligned}$$

Let us focus on the inner maximization problem, which is a convex optimization problem. So, by taking derivative with respect to $y$ and then setting that equal to zero, we can compute the maximizer as follows: 

$$
\begin{aligned}
 & \nabla_{y}\left[\sum_{i\in I}\alpha_{i}\left\{ f_{i}-\frac{1}{2L}\|y-g_{i}\|^{2}+\left\langle x-x_{i}\mid y\right\rangle \right\} \right]=0\\
\Leftrightarrow & \sum_{i\in I}\alpha_{i}\left\{ -\frac{1}{\cancel{2}}\frac{1}{L}\times\cancel{2}\times(y-g_{i})+x-x_{i}\right\} =\sum_{i\in I}\alpha_{i}\left(-\frac{1}{L}y+\frac{1}{L}g_{i}+x-x_{i}\right)=0\\
\Leftrightarrow & \sum_{i\in I}\alpha_{i}\left(\frac{1}{L}g_{i}+x-x_{i}\right)=\sum_{i\in I}\alpha_{i}\left(\frac{1}{L}y\right)=\frac{1}{L}y\underbrace{\sum_{i\in I}\alpha_{i}}_{=1}=\frac{1}{L}y\\
\Leftrightarrow & \sum_{i\in I}\alpha_{i}\left(\frac{1}{L}g_{i}+x-x_{i}\right)=\frac{1}{L}y\\
\Leftrightarrow & y^{\star}=L\Big(\sum_{i\in I}\alpha_{i}\{(x-x_{i})+\frac{1}{L}g_{i}\}\Big)\\
\therefore y^{\star} & =L\Big(x-\sum_{i\in I}\alpha_{i}\big(x_{i}-\frac{1}{L}g_{i}\big)\Big)\quad\quad\quad(4)
\end{aligned}
$$


Putting the value of $y$​​​​​ in (4) in the inner maximization problem, we have: 

$$
\begin{aligned}
 & \sup_{y\in\mathbf{R}^{d}}\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\underbrace{\|y-g_{i}\|^{2}}_{=\|y\|^{2}+\|g_{i}\|^{2}-2\langle y\mid g_{i}\rangle}+\left\langle x-x_{i}\mid y\right\rangle \Big)\\
 & =\sup_{y\in\mathbf{R}^{d}}\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\{\|y\|^{2}+\|g_{i}\|^{2}-2\langle y\mid g_{i}\rangle\}+\left\langle x-x_{i}\mid y\right\rangle \Big)\\
 & =\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\Big)+\sup_{y\in\mathbf{R}^{d}}\sum_{i\in I}\alpha_{i}\Big(-\frac{1}{2L}\|y\|^{2}+\langle y\mid\frac{1}{L}g_{i}\rangle+\left\langle x-x_{i}\mid y\right\rangle \Big)\\
 & =\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\Big)+\sup_{y\in\mathbf{R}^{d}}\sum_{i\in I}\alpha_{i}\Big(-\frac{1}{2L}\|y\|^{2}+\left\langle x-x_{i}+\frac{1}{L}g_{i}\mid y\right\rangle \Big)\\
 & =\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\Big)+\sum_{i\in I}\alpha_{i}\Big(-\frac{1}{2L}\|y^{\star}\|^{2}+\left\langle x-x_{i}+\frac{1}{L}g_{i}\mid y^{\star}\right\rangle \Big)\\
 & =\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\Big)-\frac{1}{2L}\|y^{\star}\|^{2}+\sum_{i\in I}\alpha_{i}\left\langle x-x_{i}+\frac{1}{L}g_{i}\mid y^{\star}\right\rangle \\
 & =\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\Big)-\frac{1}{2L}\|L\Big(x-\sum_{i\in I}\alpha_{i}\big(x_{i}-\frac{1}{L}g_{i}\big)\Big)\|^{2}\\
 & +\left\langle \Big(x-\sum_{i\in I}\alpha_{i}\big(x_{i}-\frac{1}{L}g_{i}\big)\Big)\mid L\Big(x-\sum_{i\in I}\alpha_{i}\big(x_{i}-\frac{1}{L}g_{i}\big)\Big)\right\rangle \\
 & =\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\Big)-\frac{L}{2}\|x-\sum_{i\in I}\alpha_{i}\big(x_{i}-\frac{1}{L}g_{i}\big)\Big)\|^{2}+L\|x-\sum_{i\in I}\alpha_{i}\big(x_{i}-\frac{1}{L}g_{i}\big)\Big)\|^{2}\\
 & =\sum_{i\in I}\alpha_{i}\Big(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\Big)+\frac{L}{2}\|x-\sum_{i\in I}\alpha_{i}\big(x_{i}-\frac{1}{L}g_{i}\big)\Big)\|^{2}\\
 & =\frac{L}{2}\|x-\sum_{i\in I}\alpha_{i}\left(x_{i}-\frac{1}{L}g_{i}\right)\|^{2}+\sum_{i\in I}\alpha_{i}\left(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\right).\quad\quad\quad(5)
\end{aligned}
$$
Using $(5)$​​​​​ in (3), we arrive at the desired interpolating function $f\in\mathcal{F}_{0,L}(\mathbf{R}^{d}):$​​​​​ $$f(x)=\min_{\alpha\in\Delta}\left[\frac{L}{2}\|x-\sum_{i\in I}\alpha_{i}\left(x_{i}-\frac{1}{L}g_{i}\right)\|^{2}+\sum_{i\in I}\alpha_{i}\left(f_{i}-\frac{1}{2L}\|g_{i}\|^{2}\right)\right],$$​​​​​ thus completing the proof.

### References.

`[1]` Taylor, Adrien B. Convex interpolation and performance estimation of first-order methods for convex optimization. Diss. Catholic University of Louvain, Louvain-la-Neuve, Belgium, 2017.

`[2]` R Tyrell Rockafellar. Convex Analysis. Princeton University Press, 1996.

`[3]` Drori, Yoel. The exact information-based complexity of smooth convex minimization. Journal of Complexity 39 (2017): 1-16.
