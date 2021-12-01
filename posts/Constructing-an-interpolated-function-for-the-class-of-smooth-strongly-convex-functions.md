\enabletheorems

@def title = "Constructing an interpolated function for the class of smooth strongly convex functions"
@def published ="December 1, 2021"
@def tags =["optimization", "lower-bounds-in-optimization"]

# Constructing an interpolated function for the class of smooth strongly convex functions

**Shuvomoy Das Gupta**

*December 1, 2021*

In this blog, we study constructing an interpolated smooth and strongly convex function from a set of points due to Yoel Drori and Adrien Taylor from`[4]`.

---
**Table of contents**
\toc

---

## Notation and notions.

All norms are 2-norm in this blog. A function $f:\mathbf{R}^{d}\to\mathbf{R}$ is $L$-smooth convex if and only $$\left(\forall x,y\in\mathbf{R}^{d}\right)\quad f(y)\geq f(x)+\langle\nabla f(x)\mid y-x\rangle+\frac{1}{2L}\|\nabla f(x)-\nabla f(y)\|^{2}.$$

On the other hand, a function is $\mu$​-strongly convex if and only if $$\left(\forall x,y\in\mathbf{R}^{d}\right)\quad f(y)\geq f(x)+\langle\nabla f(x)\mid y-x\rangle+\frac{\mu}{2}\|x-y\|^{2}\qquad (\text{SCVX})$$​ where $f^{\prime}(\cdot)$​ denotes a subgradient of $f$​ at $(\cdot)$​.

On $\mathbf{R}^{d},$ the set of all $L$-smooth convex functions is denoted by $\mathcal{F}_{0,L}(\mathbf{R}^{d})$, the set of all $\mu$-strongly convex functions is denoted by $\mathcal{F}_{\mu,\infty}(\mathbf{R}^{d})$, and the set of all $\mu$-strongly convex and $L$-smooth functions is denoted by $\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$. Finally the set of all lower-semicontinuous, proper, and convex functions is denoted by $\mathcal{F}_{0,\infty}(\mathbf{R}^{d})$.

## Helper results.

### Alternative characterization of smooth convex functions

An alternative characterization of $f\in\mathcal{F}_{0,L}(\mathbf{R}^{d})$, that we will use multiple times, is as follows.

\begin{theorem}{Alternative characterization of smooth convex functions`[1, Definition 2.6, Theorem 2.27, Theorem 2.28]`}{thm-smth-cvx-alt}

For a function $f:\mathbf{R}^{d}\to\mathbf{R}$​ the following statements are equivalent:

(a) $f\in\mathcal{F}_{0,L}(\mathbf{R}^{d})$​.

(b)$f$​ is convex and it satisfies $\left(\forall x,y\in\mathbf{R}^{d}\right)\quad f(y)\leq f(x)+\langle\nabla f(x)\mid y-x\rangle+\frac{L}{2}\|x-y\|^{2}.$​

(c) $f$ is convex and satisfies $\left(\forall x,y\in\mathbf{R}^{d}\right)\quad\|\nabla f(x)-\nabla f(y)\|\leq L\|x-y\|$.

(d) $\frac{L}{2}\|\cdot\|^{2}-f\in\mathcal{F}_{0,\infty}(\mathbf{R}^{d})$​​.
\end{theorem}

Note that, in $(b)$ and $(c)$ saying that $f$ is convex is necessary, as the inequalities only in $(b),(c)$ are not sufficient to establish $L$-smooth convexity.

We now record the notion of an interpolable function.

\begin{definition}{Interpolable function.}{def-intpl-fn}
Suppose we are given the set of triplets $\{(x_{i},g_{i},f_{i})\}_{i\in I}\subseteq\mathbf{R}^{d}\times\mathbf{R}^{d}\times\mathbf{R}$ where $I$ is a finite index set. Let $\mathcal{F}(\mathbf{R}^d)$ be some class of functions. Then the set $\{(x_{i},g_{i},f_{i} ) \}_{i\in I}$ is $\mathcal{F}(\mathbf{R}^{d})$-interpolable if and only if there exists a function $f\in\mathcal{F}(\mathbf{R}^{d})$ such that for all $i\in I$ we have $f_{i}=f(x_{i})$ and $g_{i} \in \partial f(x_{i})$.
\end{definition}

### Equivalence between smooth strongly convex and smooth convex function class

To establish our main result, we will use the following equivalence result.

\begin{lemma}{Equivalence between function classes}{lem-fun-class-equiv}

We have $f\in\mathcal{F}_{\mu,L}(\mathbf{R}^{d})\Leftrightarrow\frac{L}{2}\|\cdot\|^{2}-f\in\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})$.

\end{lemma}

**Proof.** 

\begin{dropdown}{Click to show proof}

First, we prove $f\in\mathcal{F}_{\mu,L}(\mathbf{R}^{d})\Rightarrow\frac{L}{2}\|\cdot\|^{2}-f\in\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})$​​. Define $h\coloneqq\frac{L}{2}\|\cdot\|^{2}-f$​​, so $f=\frac{L}{2}\|\cdot\|^{2}-h$​​ and $\nabla f(x)=Lx-\nabla h(x)$​​. Clearly, $f\in\mathcal{F}_{0,L}(\mathbf{R}^{d})$​​ also, so due to \theoremref{thm-smth-cvx-alt}$(a),(d),$​​we have $h\in\mathcal{F}_{0,\infty}(\mathbf{R}^{d})$​​. Because $f$​​ is $\mu$​​-strongly convex and differentiable, for any $x,y\in\mathbf{R}^{d}$​​ we have from $(\text{SCVX})$:

 $$\begin{aligned}
 & f(y)\geq f(x)+\left\langle \nabla f(x)\mid y-x\right\rangle +\frac{\mu}{2}\|x-y\|^{2}\\
\Leftrightarrow & \frac{L}{2}\|y\|^{2}-h(y)\geq\frac{L}{2}\|x\|^{2}-h(x)+\left\langle Lx-\nabla h(x)\mid y-x\right\rangle +\frac{\mu}{2}\|x-y\|^{2}\\
\Leftrightarrow & -\frac{L}{2}\|y\|^{2}+h(y)\leq-\frac{L}{2}\|x\|^{2}+h(x)-\left\langle Lx-\nabla h(x)\mid y-x\right\rangle -\frac{\mu}{2}\|x-y\|^{2}\\
\Leftrightarrow & h(y)\leq\frac{L}{2}\|y\|^{2}-\frac{L}{2}\|x\|^{2}+h(x)-L\left\langle x\mid y-x\right\rangle +\left\langle \nabla h(x)\mid y-x\right\rangle -\frac{\mu}{2}\|x-y\|^{2}\\
 & \quad\quad\quad=\frac{L}{2}\|y\|^{2}-\frac{L}{2}\|x\|^{2}+h(x)-L\left\langle x\mid y\right\rangle +L\|x\|^{2}+\left\langle \nabla h(x)\mid y-x\right\rangle -\frac{\mu}{2}\|x-y\|^{2}\\
 & \quad\quad\quad=\frac{L}{2}\left(\|y\|^{2}+\|x\|^{2}-2\left\langle x\mid y\right\rangle \right)+h(x)+\left\langle \nabla h(x)\mid y-x\right\rangle -\frac{\mu}{2}\|x-y\|^{2}\\
 & \quad\quad\quad=h(x)+\left\langle \nabla h(x)\mid y-x\right\rangle +\frac{L}{2}\|x-y\|^{2}-\frac{\mu}{2}\|x-y\|^{2}\\
 & \quad\quad\quad=h(x)+\left\langle \nabla h(x)\mid y-x\right\rangle +\frac{L-\mu}{2}\|x-y\|^{2}.\end{aligned}$$

So we have proven that $h$ is convex and for any $x,y\in\mathbf{R}^{d},$ $$h(y)\leq h(x)+\left\langle \nabla h(x)\mid y-x\right\rangle +\frac{L-\mu}{2}\|x-y\|^{2}$$ which this is equivalent to saying that $h\in\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})$ due to \theoremref{thm-smth-cvx-alt}$(a),(b)$.

Next, we prove $h\coloneqq\frac{L}{2}\|\cdot\|^{2}-f\in\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})\Rightarrow f\in\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$. Due to \theoremref{thm-smth-cvx-alt}$(a),(b)$, we have: $h$ convex and for any $x,y\in\mathbf{R}^{d},$ 

$$\begin{aligned}
 & h(y)\leq h(x)+\left\langle \nabla h(x)\mid y-x\right\rangle +\frac{L-\mu}{2}\|x-y\|^{2}\\
\Leftrightarrow & \frac{L}{2}\|y\|^{2}-f(y)\leq\frac{L}{2}\|x\|^{2}-f(x)+\left\langle Lx-\nabla f(x)\mid y-x\right\rangle +\frac{L-\mu}{2}\|x-y\|^{2}\\
\Leftrightarrow & -\frac{L}{2}\|y\|^{2}+f(y)\geq-\frac{L}{2}\|x\|^{2}+f(x)-\left\langle Lx-\nabla f(x)\mid y-x\right\rangle -\frac{L-\mu}{2}\|x-y\|^{2}\\
\Leftrightarrow & f(y)\geq f(x)+\left\langle \nabla f(x)\mid y-x\right\rangle +\frac{L}{2}\|y\|^{2}-\frac{L}{2}\|x\|^{2}-L\left\langle x\mid y\right\rangle +L\|x\|^{2}-\frac{L-\mu}{2}\|x-y\|^{2}\\
 & \quad\quad\quad=f(x)+\left\langle \nabla f(x)\mid y-x\right\rangle +\frac{L}{2}\underbrace{\left(\|y\|^{2}+\|x\|^{2}-2\left\langle x\mid y\right\rangle \right)}_{=\|x-y\|^{2}}-\frac{L-\mu}{2}\|x-y\|^{2}\\
 & \quad\quad\quad=f(x)+\left\langle \nabla f(x)\mid y-x\right\rangle +\frac{\mu}{2}\|x-y\|^{2},\end{aligned}$$

*i.e.*, we have shown that for all $x,y\in\mathbf{R}^{d}$​​​ $$f(y)\geq f(x)+\left\langle \nabla f(x)\mid y-x\right\rangle +\frac{\mu}{2}\|x-y\|^{2},$$​​​ hence by defintion $f$​​​ is $\mu$​​​-strongly convex. Finally, we have $$\begin{aligned}
\|\nabla f(x)-\nabla f(y)\| & =\|Lx-\nabla h(x)-Ly+\nabla h(y)\|\\
 & =\|L(x-y)+(\nabla h(y)-\nabla h(x))\|\\
 & \overset{a)}{\leq}L\|x-y\|+\|\nabla h(y)-\nabla h(x)\|\\
 & \overset{b)}{\leq}L\|x-y\|+(L-\mu)\|x-y\|\\
 & =L\|x-y\|,\end{aligned}$$

where $a)$ follows from triangle inequallity and $b)$ follows from \theoremref{thm-smth-cvx-alt}$(a),(c)$. So, $f$ is $L$-smooth besides being $\mu$-strongly convex, *i.e.*, $f\in\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$. This completes the proof. 

\end{dropdown}
■

### Interpolation equivalence for smooth strongly convex and smooth convex functions

Also, we are going to use the following interpolation result from`[1, Theorem 3.8]`.

\begin{theorem}{Interpolation equivalence for smooth strongly convex and smooth convex functions}{thm-int-eqvl}

Suppose we are given the set of triplets set of triplets $\{(x_{i},g_{i},f_{i})\}_{i\in I}\subseteq\mathbf{R}^{d}\times\mathbf{R}^{d}\times\mathbf{R}$ where $I$​ is a finite index set. Then the following are equivalent.

(a) $\{(x_{i},g_{i},f_{i})\}_{i\in I}$​ is $\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$​​-interpolable.

(b) $\{(\widetilde{x}_{i},\widetilde{g}_{i},\widetilde{f}_{i})\}_{i\in I}\coloneqq\{(x_{i},Lx_{i}-g_{i},\frac{L}{2}\|x_{i}\|^{2}-f_{i})\}_{i\in I}$​ is $\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})$​​-interpolable.
\end{theorem}

**Proof.** 

\begin{dropdown}{Click to show proof}

First, we prove $(a)\Rightarrow(b).$​ If $(a)$​ holds, then by definition there exists a function $f\in\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$​ that satisfies for all $i\in I$: $f(x_{i})=f_{i}$​ and $\nabla f(x_{i})=g_{i}$​. If we define, $h=(L/2)\|\cdot\|^{2}-f,$​ then $h\in\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})$​ due to \lemmaref{lem-fun-class-equiv}, and for all $i\in I$​ we have $h(x_{i})=\frac{L}{2}\|x_{i}\|^{2}-f(x_{i})=\frac{L}{2}\|x_{i}\|^{2}-f_{i}$​ and $\nabla h(x_{i})=Lx_{i}-\nabla f(x_{i})=Lx_{i}-g_{i}$​. This proves $(b)$​.

Next, we prove $(b)\Rightarrow(a)$. If $(b)$ holds, there exists a function $h\in\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$ that satisfies for all $i\in I$: $h(x_{i})=\frac{L}{2}\|x_{i}\|^{2}-f_{i}$ and $\nabla h(x_{i})=Lx_{i}-g_{i}$. If we define, $f=(L/2)\|\cdot\|^{2}-h$, then $f\in\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})$ due to  \lemmaref{lem-fun-class-equiv}, and for all $i\in I$ we have $f(x_{i})=\frac{L}{2}\|x_{i}\|^{2}-h(x_{i})=\frac{L}{2}\|x_{i}\|^{2}-(\frac{L}{2}\|x_{i}\|^{2}-f_{i})=f_{i}$ and $\nabla f(x_{i})=Lx_{i}-\nabla h(x_{i})=Lx_{i}-(Lx_{i}-g_{i})=g_{i}.$​ This proves (a). 

\end{dropdown}
■


### Drori's interpolated function for smooth convex function class
Finally, we present the following interpolation result due to Yoel Drori from`[3]`. I showed a detailed proof of this result in the previous blog post: [here](https://shuvomoy.github.io/blogs/posts/Constructing-an-interpolated-function-for-the-class-of-smooth-convex-functions/).

\begin{theorem}{Interpolation of smooth convex functions.}{Thm-interpolation-of-smooth}
If $\{(\widetilde{x}_{i},\widetilde{g}_{i},\widetilde{f}_{i})\}_{i\in I}\subseteq\mathbf{R}^{d}\times\mathbf{R}^{d}\times\mathbf{R}$​ is $\mathcal{F}_{0,L}(\mathbf{R}^{d})$​-interpolable with $L>0$​, then one interpolated function $\widetilde{f}\in\mathcal{F}_{0,L}(\mathbf{R}^{d})$​ that interpolates $\{(\widetilde{x}_{i},\widetilde{g}_{i},\widetilde{f}_{i})\}_{i\in I}$​ is: 
$$\widetilde{f}(y)=\min_{\alpha\in\Delta}\left[\frac{L}{2}\|y-\sum_{i\in I}\alpha_{i}\left(\widetilde{x}_{i}-\frac{1}{L}\widetilde{g}_{i}\right)\|^{2}+\sum_{i\in I}\alpha_{i}\left(\widetilde{f}_{i}-\frac{1}{2L}\|\widetilde{g}_{i}\|^{2}\right)\right],$$​ 
where $\Delta=\{\beta\in\mathbf{R}^{\vert I\vert}\mid\beta\geq0,\sum_{i\in I}\beta_{i}=1\}.$​
\end{theorem}

## Main result.

Now we are in a position to state our main result due to Yoel Drori and Adrien Taylor from`[4, Theorem 1]` followed by its proof. Their proof is direct and the proof we present here is based on interpolation argument.



\begin{theorem}{Interpolation of smooth strongly convex functions}{Thm-interpolation-of-smth-scvx}

If $\{(x_{i},g_{i},f_{i})\}_{i\in I}\subseteq\mathbf{R}^{d}\times\mathbf{R}^{d}\times\mathbf{R}$​ is $\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$​-interpolable with $L>0$​, then one interpolation function $f\in\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$​ that interpolates $\{(x_{i},g_{i},f_{i})\}_{i\in I}$​ is: 
$$\begin{align*}
f(y) & :=\max_{\alpha\in\Delta}\Big[\frac{L}{2}\|y\|^{2}-\frac{L-\mu}{2}\|y-\frac{1}{L-\mu}\sum_{i}\alpha_{i}(g_{i}-\mu x_{i})\|^{2}\\
 & \quad  +\sum_{i\in I}\alpha_{i}(f_{i}+\frac{1}{2(L-\mu)}\|g_{i}-Lx_{i}\|^{2}-\frac{L}{2}\|x_{i}\|^{2})\Big],
\end{align*}$$
where $\Delta=\{\beta\in\mathbf{R}^{\vert I\vert}\mid\beta\geq0,\sum_{i\in I}\beta_{i}=1\}.$

\end{theorem}

**Proof.** The proof sketch is as follows comprising two steps. 

1.  First, we find an interpolated function $\widetilde{f}\in\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})$ that interpolates $\{(\widetilde{x}_{i},\widetilde{g}_{i},\widetilde{f}_{i})\}_{i\in I}\coloneqq\{(x_{i},Lx_{i}-g_{i},\frac{L}{2}\|x_{i}\|^{2}-f_{i})\}_{i\in I}$ using  \theoremref{Thm-interpolation-of-smooth}.

2.  Then due to the proof of \theoremref{thm-int-eqvl}, the function $f=(L/2)\|\cdot\|^{2}-\widetilde{f}$ will be in $\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$ and will interpolate $\{(x_{i},g_{i},f_{i})\}_{i\in I}$.

\(1\) From \theoremref{Thm-interpolation-of-smooth} recall that, if $\{(\widetilde{x}_{i},\widetilde{g}_{i},\widetilde{f}_{i})\}_{i\in I}\subseteq\mathbf{R}^{d}\times\mathbf{R}^{d}\times\mathbf{R}$​ is $\mathcal{F}_{0,L-\mu}(\mathbf{R}^{d})$​-interpolable, then one interpolation function $\widetilde{f}\in\mathcal{F}_{0,L}(\mathbf{R}^{d})$​ that interpolates $\{(\widetilde{x}_{i},\widetilde{g}_{i},\widetilde{f}_{i})\}_{i\in I}$​ is: 

$$\widetilde{f}(y)=\min_{\alpha\in\Delta}\left[\frac{L-\mu}{2}\|y-\sum_{i\in I}\alpha_{i}\left(\widetilde{x}_{i}-\frac{1}{L-\mu}\widetilde{g}_{i}\right)\|^{2}+\sum_{i\in I}\alpha_{i}\left(\widetilde{f}_{i}-\frac{1}{2(L-\mu)}\|\widetilde{g}_{i}\|^{2}\right)\right],$$

where $\Delta=\{\beta\in\mathbf{R}^{\vert I\vert}\mid\beta\geq0,\sum_{i\in I}\beta_{i}=1\}.$​​ In our setup, we have $\{(\widetilde{x}_{i},\widetilde{g}_{i},\widetilde{f}_{i})\}_{i\in I}\coloneqq\{(x_{i},Lx_{i}-g_{i},\frac{L}{2}\|x_{i}\|^{2}-f_{i})\}_{i\in I}$​​, so lets put that in the last equation:

$$\begin{aligned}
\widetilde{f}(y) & =\min_{\alpha\in\Delta}\left[\frac{L-\mu}{2}\|y-\sum_{i\in I}\alpha_{i}\Big(\underbrace{x_{i}-\frac{1}{L-\mu}(Lx_{i}-g_{i})}_{=x_{i}(1-\frac{L}{L-\mu})+\frac{1}{L-\mu}g_{i}=-\frac{\mu}{L-\mu}x_{i}+\frac{1}{L-\mu}g_{i}=\frac{1}{L-\mu}(g_{i}-\mu x_{i})}\Big)\|^{2}+\sum_{i\in I}\alpha_{i}\left(\frac{L}{2}\|x_{i}\|^{2}-f_{i}-\frac{1}{2(L-\mu)}\|Lx_{i}-g_{i}\|^{2}\right)\right]\\
 & =\min_{\alpha\in\Delta}\left[\frac{L-\mu}{2}\|y-\frac{1}{L-\mu}\sum_{i\in I}\alpha_{i}(g_{i}-\mu x_{i})\|^{2}-\sum_{i\in I}\alpha_{i}\left(f_{i}+\frac{1}{2(L-\mu)}\|g_{i}-Lx_{i}\|^{2}-\frac{L}{2}\|x_{i}\|^{2}\right)\right].\end{aligned}$$

(2) Hence, $f=(L/2)\|\cdot\|^{2}-\widetilde{f}$​​​, which will be in $\mathcal{F}_{\mu,L}(\mathbf{R}^{d})$​​​ and will interpolate $\{(x_{i},g_{i},f_{i})\}_{i\in I}$​​​ has the following form: 

$$\begin{aligned}
f(y) & =(L/2)\|y\|^{2}-\widetilde{f}(y)\\
 & =(L/2)\|y\|^{2}-\min_{\alpha\in\Delta}\left[\frac{L-\mu}{2}\|y-\frac{1}{L-\mu}\sum_{i\in I}\alpha_{i}(g_{i}-\mu x_{i})\|^{2}-\sum_{i\in I}\alpha_{i}\left(f_{i}+\frac{1}{2(L-\mu)}\|g_{i}-Lx_{i}\|^{2}-\frac{L}{2}\|x_{i}\|^{2}\right)\right]\\
 & \overset{a)}{=}(L/2)\|y\|^{2}+\max_{\alpha\in\Delta}\left[-\frac{L-\mu}{2}\|y-\frac{1}{L-\mu}\sum_{i\in I}\alpha_{i}(g_{i}-\mu x_{i})\|^{2}+\sum_{i\in I}\alpha_{i}\left(f_{i}+\frac{1}{2(L-\mu)}\|g_{i}-Lx_{i}\|^{2}-\frac{L}{2}\|x_{i}\|^{2}\right)\right]\\
 & =\max_{\alpha\in\Delta}\left[(L/2)\|y\|^{2}-\frac{L-\mu}{2}\|y-\frac{1}{L-\mu}\sum_{i\in I}\alpha_{i}(g_{i}-\mu x_{i})\|^{2}+\sum_{i\in I}\alpha_{i}\left(f_{i}+\frac{1}{2(L-\mu)}\|g_{i}-Lx_{i}\|^{2}-\frac{L}{2}\|x_{i}\|^{2}\right)\right],\end{aligned}$$

where $a)$​​​ uses $\min(\cdot)=-\max(-\cdot)$​​​. This completes the proof. ■

## References.

`[1]` Taylor, Adrien B. Convex interpolation and performance estimation of first-order methods for convex optimization. Diss. Catholic University of Louvain, Louvain-la-Neuve, Belgium, 2017.

`[2]` R Tyrell Rockafellar. Convex Analysis. Princeton University Press, 1996.

`[3]` Drori, Yoel. The exact information-based complexity of smooth convex minimization. Journal of Complexity 39 (2017): 1-16.

`[4]` Drori, Yoel and Taylor, Adrien, 2021. On the oracle complexity of smooth strongly convex minimization. arXiv preprint arXiv:2101.09740.
