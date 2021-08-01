@def title = "Geometry of a prox-regular set"
@def published ="July 30, 2021"
@def tags =["optimization"]

# Geometry of a prox-regular set

**Shuvomoy Das Gupta**

*July 30, 2021*

In this blog, we provide a proof about certain geometric properties of a prox-regular set based on Proposition 3.1 of the paper "Local Differentiability of Distance Functions" by R. A. Poliquin, R. T. Rockafellar and L. Thibault (link to pdf [here](https://www.ams.org/journals/tran/2000-352-11/S0002-9947-00-02550-2/S0002-9947-00-02550-2.pdf)). 

---

\toc

---

We consider nonempty, closed set $\mathcal{X}$ in a finite-dimensional vector space, which is *prox-regular* at a point $\bar{x}\in\mathcal{X}$, *i.e.*, Euclidean projection onto the set from $\bar{x}$ is single-valued in some neighborhood of that point. Throughout the blog, the underlying space is a finite-dimensional vector space over the reals. To state the result and to prove it, we introduce certain notation and notions.

## Notation and notions

### Notation.

An open ball with center $c$ and radius $r$ is denoted by 
$$
B(c;r)=\left\{ x\mid\|x-c\| < r\right\},
$$
where the norm is is the 2 norm throughout the blog. 

#### Indicator function.

The indicator function of $\mathcal{X}$ is denoted by $\iota,$ with 

$$\iota(x)=\begin{cases}
0, & x\in\mathcal{X}\\
\infty, & x\notin\mathcal{X}.
\end{cases}$$

#### Projection.

A Euclidean projection from a point $x$ onto this set $\mathcal{X}$ is denoted by $\Pi(x)$, and one arbitrary projection in $\Pi(x)$ is denoted by $\pi_{x}$. It is defined as: 
$$\pi_{x}\in\Pi(x)=\mathop{\textrm{argmin}_{y\in\mathcal{X}}\|x-y\|^{2}.}\quad(\textrm{Projx})$$​ 
The distance function at the point $x$ is denoted by $d(x)=\|x-\pi_{x}\|$, where the norm is the 2 norm.

### The distance squared function is locally Lipschitz.

A function is locally Lipschitz, if at any point, there is a neighborhood around that point where the function is Lipschitz continuous.

For any closed, nonempty set its distance function is 1-Lipschitz continuous everywhere \citep{Rockafellar09}[Example 9.6]. This implies that for any $\beta$, $d^{2}+\beta\|x\|^{2}$ is locally Lipschitz on any point $x$. To see that, consider $x\in\mathcal{B}$ where $\mathcal{B}$ is some closed ball around $x$. Then, for any $x,y\in\mathcal{B},$ 
$$
\begin{align*}
 & |d^{2}(x)+\beta\|x\|^{2}-d^{2}(y)-\beta\|y\|^{2}|\\
\leq & |\left(d(x)-d(y)\right)\left(d(x)+d(y)\right)|+|\beta||\underbrace{\|x\|^{2}-\|y\|^{2}}_{=(\|x\|+\|y\|)(\|x\|-\|y\|)}|\\
\leq & \underbrace{|d(x)+d(y)|}_{\leq M_{1}}\underbrace{|d(x)-d(y)|}_{\leq\|x-y\|}\\
 & +|\beta|\underbrace{|\|x\|+\|y\||}_{\leq M_{2}}\underbrace{|\|x\|-\|y\||}_{\leq\|x-y\|}\\
\leq & (M_{1}+M_{2}\beta)\|x-y\|,\quad(\textrm{DistSqdLcLip})
\end{align*}
$$
where in the second line we have used (i) the existence of constant $M_{1},M_{2}>0$​ due to the points $x,y\in\mathcal{B}$, and (ii) the reverse triangle inequality.

### Proximal normal cone.

The proximal normal cone of $\mathcal{X}$ at a point $x\in\mathcal{X}$ is defined as follows  \citep{Rockafellar09}[Example 6.16]: 
$$
v\in N(x)\Leftrightarrow\exists_{\tau>0}\;x\in\Pi(x+\tau v),\quad(\textrm{PrxNrmlCn})
$$
which also implies that
$$
\forall_{\widetilde{\tau}\in(0,\tau)}\;\Pi(x+\widetilde{\tau}v)=\{x\}.
$$
### Frechet subdifferential and Clarke subdifferential. 

For any function $f$ (not necessarily convex), its Frechet subdifferential $\partial f$ at a point $x$ is defined as follows \citep{Correa92}[Definition 2.5]: 
$$v\in\partial f(x)\Leftrightarrow\liminf_{y\to0}\frac{f(x+y)-f(x)-\left\langle v\mid y\right\rangle }{\|y\|}\geq0.$$​

On the other hand, The Clarke subdifferential of a locally Lipschitz function $f$​ is defined as follows: 
$$
u\in\partial^{\textrm{Clarke}}f(x)\Leftrightarrow\forall_{d}\;\left[\limsup_{y\to x,t\downarrow0}\frac{f(y+td)-f(y)}{t}\right]\geq\left\langle u\mid d\right\rangle .
$$

For a locally Lipschitz function the Clarke subdifferential is nonempty everywhere \citep{Correa92}[Property 2.2]. 

#### Finding Frechet subgradient of an infimal convolution function.

Define the infimal convolution function between two lower-semicontinuous functions $f$ and $g$ as follows \citep{Correa92}[Lemma 3.6]:
$$\begin{aligned}
(f\square g)(x) & =\inf_{y}\left\{ f(y)+g(x-y)\right\} .\end{aligned}$$ 
Also, define 
$$y_{x}^{\star}\in\mathop{\textrm{argmin}}_{y}\left\{ f(y)+g(x-y)\right\},$$
which can be empty. If $y_{x}^{\star}$ exists, then
$$\partial(f\square g)(x)\subseteq\partial f(y_{x}^{\star})\cap\partial g(x-y_{x}^{\star}).\quad(\textrm{FrechSubRule})$$

#### An interesting implication of (FrechSubRule).

Note that, 
$$d^{2}(x)=\min_{y\in\mathcal{X}}\|x-y\|^{2}=\min_{y}\left\{ \iota(y)+\|x-y\|^{2}\right\} =(\iota\square\|\cdot\|^{2})(x),$$​ 
where by definition: 
$$\pi_{x}=\mathop{\textrm{argmin}}_{y\in\mathcal{X}}\|x-y\|^{2}=\mathop{\textrm{argmin}}_{y}\left\{ \iota(y)\square\|x-y\|^{2}\right\},$$
Hence, recalling that $\partial\|x\|^{2}=2x$, we can find a Frechet subgradient of $d^{2}$ using (FrechSubRule):
$$\begin{aligned}
\partial d^{2}(x) & =\partial(\iota\square\|\cdot\|^{2})(x)\\
\subseteq & \partial\iota(\pi_{x})\cap\underbrace{\partial\left[\|x-\pi_{x}\|^{2}\right]}_{=2(x-\pi_{x})}\\
 & =\partial\iota(\pi_{x})\cap\{2(x-\pi_{x})\},\end{aligned}$$
hence if $\partial d^{2}(x)\neq\emptyset$, then $\partial d^{2}(x)=\{2(x-\pi_{x})\}$. In other words,
$$\forall_{x:\partial d^{2}(x)\neq\emptyset}\quad\partial d^{2}(x)=2(x-\pi_{x}).\quad(\textrm{FrechSubDistSqd})$$

#### Proving local convexity via Frechet subdifferential.

If a locally Lipschitz function $f$​ has its Frechet subdifferential $\partial f$​ monotone on $\{(x,u) \in \mathbf{gra}\partial f \mid x\in A\}$​, then $f$​ is convex on $A$​ \citep{Correa92}[Theorem 3.8, Remark after Property 2.7]. In other words, for a locally Lipschitz function $f$​
$$\forall_{(x,u),(y,v)\in\mathbf{gra}\partial f:x,y\in A}\;\left\langle u-v\mid x-y\right\rangle \geq0\Rightarrow f:\textrm{convex on }A\quad(\textrm{LocCvx}),$$​
where $\mathbf{gra}\partial f=\left\{ (x,u)\mid u\in\partial f(x)\right\}.$

Furthermore, if $f$ is locally Lipschitz, then monotonicity of $\partial f$ is equivalent to the monotonicity of $\partial^{\textrm{Clarke}}f$, so proving either is fine to establish convexity \citep{Correa92}[Remark after Property 2.7]. 

### Cocoercive operator.

An operator $A$ is $\beta$-cocoercive on $S$ if
$$\forall_{(x,u),(y,v)\in S\cap\mathbf{gra}A}\quad\left\langle x-y\mid u-v\right\rangle \geq\beta\|u-v\|^{2}.\quad(\textrm{CocrcvOpt})$$
A $\beta$-cocoercive operator on $S$ is also $\frac{1}{\beta}$-Lipschitz on $S$.

### Locally hypomonotone operator.

An operator $A$ is $\rho$-hypomonotone on a set $S$ if
$$\forall_{(x,u),(u,v)\in S\cap\mathbf{gra}A}\quad\left\langle u-v\mid x-y\right\rangle +\rho\|x-y\|^{2}\geq0.$$

#### Proximal normal cone of a prox-regular set is locally hypomonotone.

The following result is a restatement of Equation (3.1) of \citep{Poliquin00}. If $\mathcal{X}$ is prox-regular at $\bar{x}\in\mathcal{X}$, then there exists some $\rho>0$ and $R>0$ such that the proximal normal cone is a $\rho$-hypomonotone operator on an $R$-neighborhood of $(\bar{x},0)$, defined by
$$\mathcal{V}_{R}(\bar{x},0)=\left\{ (p,u)\mid\|p-\bar{x}\|\leq R,\|u-0\|\leq R\right\},\quad \textrm{(VR)}$$
*i.e.*, 
$$\begin{aligned}
\forall_{(p,u),(q,v)\in\mathcal{V}_{R}(\bar{x},0)\cap\mathbf{gra}N} & \quad\left\langle u-v\mid p-q\right\rangle +\rho\|p-q\|^{2}\geq 0. \quad(\textrm{HypMntn})\end{aligned}$$
Now we are in a position to state the main result regarding geometry of a prox-regular set and prove it.

## Geoemetry of a prox-regular set.


If $\mathcal{X}$ is set that is prox-regular at $\bar{x}\in\mathcal{X}$, then there exist some $\rho>0,R>0$ such that for any $\lambda$ satisfying $\lambda\in(0,2)$ and $\lambda\leq\rho$ we have the following properties :

\(i\) for any $x\in B(\bar{x};\frac{\lambda R}{2\rho})$​, we have $\left(\pi_{x},\frac{2\rho}{\lambda}(x-\pi_{x})\right)\in\mathcal{V}_{R}(\bar{x},0)\cap\mathbf{gra}N$​, where $\mathcal{V}_R(\bar x,0)$ is the same set defined in (VR).

\(ii\) the projection operator $\Pi$ is $\frac{2}{2-\lambda}$-Lipschitz continuous on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$,

\(iii\) on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$​, the function $\phi_{\lambda}=d^{2}+\frac{\lambda}{2-\lambda}\|\cdot\|^{2}$ is convex and differentiable, with the derivative given by $\nabla\phi_{\lambda}(x)=2(x-\pi_{x})+2\frac{2}{2-\lambda}x.$​

## Proof. 


### Proof to (i).

First, we show that for for any $x\in B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$, $\|\pi_{x}-\bar{x}\|\leq R$. For $x\in B\left(\bar{x},\frac{\lambda R}{2\rho}\right),$ we have
$$\begin{aligned}
 & \|\pi_{x}-\bar{x}\|\\
= & \|\pi_{x}-x+x-\bar{x}\|\\
\overset{a)}{\leq} & \underbrace{\|\pi_{x}-x\|}_{=d(x)\leq\|\bar{x}-x\|<(\lambda R/2\rho)}+\underbrace{\|x-\bar{x}\|}_{<(\lambda R/2\rho)}\\
< & \frac{\lambda R}{\rho}\\
\leq & R,\end{aligned}$$
where in the last line we have used $\lambda\leq\rho$.

Set $\tau\coloneqq\frac{\lambda}{2\rho}>0$, then
$$\Pi\left(\pi_{x}+\tau\left\{ \frac{2\rho}{\lambda}(x-\pi_{x})\right\} \right)=\Pi\left(x\right)\ni\pi_{x},$$
where the last inclusion follows from (Projx). Thus $\frac{2\rho}{\lambda}(x-\pi_{x})\in N(\pi_{x})$ by (PrxNrmlCn). Also, from (DistEq1), 
$$\begin{aligned}
\|\frac{2\rho}{\lambda}(x-\pi_{x})\|< & R,\end{aligned}$$

thus completing proof to (i).

### Proof to (ii).

From (i), we have for any $x,y\in B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$, we have $\left(\pi_{x},\frac{2\rho}{\lambda}(x-\pi_{x})\right),\left(\pi_{y},\frac{2\rho}{\lambda}(y-\pi_{y})\right)\in\mathcal{V}_{R}(\bar{x},0)\cap\mathbf{gra}N$, and using (HypMntn), we have
$$\begin{aligned}
0 & \leq\left\langle \frac{2\rho}{\lambda}(x-\pi_{x})-\frac{2\rho}{\lambda}(y-\pi_{y})\mid\pi_{x}-\pi_{y}\right\rangle +\rho\|\pi_{x}-\pi_{y}\|^{2}\\
 & =\left\langle \frac{2\rho}{\lambda}(x-y)-\frac{2\rho}{\lambda}(\pi_{x}-\pi_{y})\mid\pi_{x}-\pi_{y}\right\rangle +\rho\|\pi_{x}-\pi_{y}\|^{2}\\
 & =\frac{2\rho}{\lambda}\left\langle x-y\mid\pi_{x}-\pi_{y}\right\rangle \underbrace{-\frac{2\rho}{\lambda}\|\pi_{x}-\pi_{y}\|^{2}+\rho\|\pi_{x}-\pi_{y}\|^{2}}_{=-\frac{2\rho}{\lambda}\left(1-\frac{\lambda}{2}\right)\|\pi_{x}-\pi_{y}\|^{2}}\\
 & =\frac{2\rho}{\lambda}\left\langle x-y\mid\pi_{x}-\pi_{y}\right\rangle -\frac{2\rho}{\lambda}\left(1-\frac{\lambda}{2}\right)\|\pi_{x}-\pi_{y}\|^{2}\\
\Leftrightarrow & \left(1-\frac{\lambda}{2}\right)\|\pi_{x}-\pi_{y}\|^{2}\leq\left\langle x-y\mid\pi_{x}-\pi_{y}\right\rangle .\end{aligned}$$
So we have shown that 
$$\forall_{(x,\pi_{x}),(y,\pi_{y})\in\mathbf{gra}\Pi\cap B\left(\bar{x},\frac{\lambda R}{2\rho}\right)}\quad\left\langle x-y\mid\pi_{x}-\pi_{y}\right\rangle \geq\left(1-\frac{\lambda}{2}\right)\|\pi_{x}-\pi_{y}\|^{2},$$ 
hence the projection operator $\Pi$ is $\left[(2-\lambda)/2\right]$-cocoercive on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$ from (CocrcvOpt). Because a $\beta$-cocoercive operator on $S$ is also $\frac{1}{\beta}$-Lipschitz on $S$, we have the projection operator $\Pi$ being $\left[2/(2-\lambda)\right]$-Lipschitz continuous on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$, *i.e.*,
$$\forall_{(x,\pi_{x}),(y,\pi_{y})\in\mathbf{gra}\Pi\cap B\left(\bar{x},\frac{\lambda R}{2\rho}\right)}\quad\|\pi_{x}-\pi_{y}\|\leq\frac{2}{2-\lambda}\|x-y\|.\quad(\textrm{LipCont})$$
This proves (ii).

### Proof to (iii).

Take a point $x\in B\left(\bar{x},\frac{\lambda R}{2\rho}\right).$Due to (ii), the projection operator $\Pi(x)$ is single-valued on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$.

From (DistSqdLcLip) $\phi_{\lambda}$ is locally Lipschitz, so we will employ (LocCvx) to prove its convexity on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$. First, we will show that on the set 
$$S=\left\{ (x,u)\mid(x,u)\in\mathbf{gra}\partial\phi_{\lambda},x\in B\left(\bar{x},\frac{\lambda R}{2\rho}\right)\right\}$$
the operator $\partial\phi_{\lambda}$​ is monotone, which will help in proving that $\phi_{\lambda}$​ is convex and differentiable on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$​.

Recall from (FrechSubDistSqd) that 
$$\forall_{x:\partial d^{2}(x)\neq\emptyset}\quad\partial d^{2}(x)=2(x-\pi_{x}).$$​
Consider two points $(x,u),(y,v)\in S$​. Without loss of generality, we can assume that $\partial d^{2}(x),\partial d^{2}(y)$​ are nonempty, because for the empty case (iii) is vacuously true. Then from (FrechSubDistSqd) we have $\partial d^{2}(x)=2(x-\pi_{x})$​ and $\partial d^{2}(y)=2(y-\pi_{y})$​. On these points, 
$$\partial\phi_{\lambda}(x)=\partial d^{2}(x)+2\frac{\lambda}{2-\lambda}(x),\textrm{ and }\partial\phi_{\lambda}(y)=\partial d^{2}(y)+2\frac{\lambda}{2-\lambda}(y).$$​
We want to show that for any such $(x,u),(y,v)\in S$​, we have
$$\begin{aligned}
0 & \leq\left\langle \partial\phi_{\lambda}(x)-\partial\phi_{\lambda}(y)\mid x-y\right\rangle \\
 & =\left\langle \partial d^{2}(x)+2\frac{\lambda}{2-\lambda}(x)-\partial d^{2}(y)-2\frac{\lambda}{2-\lambda}(y)\mid x-y\right\rangle \\
 & =\left\langle \partial d^{2}(x)-\partial d^{2}(y)\mid x-y\right\rangle +2\frac{\lambda}{2-\lambda}\|x-y\|^{2}.\quad(\textrm{GoalA})\end{aligned}$$​
To prove (GoalA), first we note that 
$$\begin{aligned}
\left\langle \partial d^{2}(x)-\partial d^{2}(y)\mid x-y\right\rangle  & =\left\langle 2(x-\pi_{x})-2(y-\pi_{y})\mid x-y\right\rangle \\
 & =2\left\langle (x-y)-(\pi_{x}-\pi_{y})\mid x-y\right\rangle \\
 & =2\|x-y\|^{2}-2\left\langle \pi_{x}-\pi_{y}\mid x-y\right\rangle .\quad(\textrm{EqPart1})\end{aligned}$$​
Now, using Cauchy--Schwarz inequality, we have 
$$\begin{aligned}
\left\langle \pi_{x}-\pi_{y}\mid x-y\right\rangle  & \leq\underbrace{\|\pi_{x}-\pi_{y}\|}_{\leq\frac{2}{2-\lambda}\|x-y\|}\|x-y\|\\
 & \leq\frac{2}{2-\lambda}\|x-y\|^{2}\\
\Rightarrow-2\left\langle \pi_{x}-\pi_{y}\mid x-y\right\rangle  & \geq-\frac{4}{2-\lambda}\|x-y\|^{2},\end{aligned}$$​ 
and putting this in (EqPart1), we have 
$$\begin{aligned}
\left\langle \partial d^{2}(x)-\partial d^{2}(y)\mid x-y\right\rangle  & \geq2\|x-y\|^{2}-\frac{4}{2-\lambda}\|x-y\|^{2}\\
 & \geq\left(2-\frac{4}{2-\lambda}\right)\|x-y\|^{2}\\
 & =\frac{-2\lambda}{2-\lambda}\|x-y\|^{2}\\
\Rightarrow\left\langle \partial d^{2}(x)-\partial d^{2}(y)\mid x-y\right\rangle +2\frac{\lambda}{2-\lambda}\|x-y\|^{2} & \geq0,\end{aligned}$$​
thus reaching (GoalA). So, we have shown that on $S$​, $\partial\phi_{\lambda}$​ is monotone on $S$​. As $\phi_{\lambda}$​ is locally Lipschitz, it means that $\partial^{\textrm{Clarke}}f$​ is monotone on $S$​, and due to (LocCvx), we have $\phi_{\lambda}$​ convex on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$​. This further implies that, for any $x$​ in $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$​, $\partial\phi_{\lambda}(x)=\partial^{\textrm{Clarke}}\phi_{\lambda}(x),$​ due to the locally Lipschitz nature of $\phi_{\lambda},$​ it has nonempty Clarke subdifferential everywhere \citep{Correa92}[Property 2.2]. So, all points in $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$​ is Frechet subdifferentiable, *i.e.,* for any $x$​ in $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$​, we have $\partial\phi_{\lambda}(x)=\partial^{\textrm{Clarke}}\phi_{\lambda}(x)\neq\emptyset$​, and as we have shown before, on those  it is in fact differentiable with the gradient 

$$\partial\phi_{\lambda}(x)=2(x-\pi_{x})+2\frac{2}{2-\lambda}x.$$​ 

Thus, on $B\left(\bar{x},\frac{\lambda R}{2\rho}\right)$​, the function $\phi_{\lambda}$​  is convex and differentiable. 

## References

* \biblabel{Rockafellar09}{Rockafellar and Wets (2009)} Rockafellar, R. Tyrrell, and Roger J-B. Wets. *Variational analysis*. Vol. 317. Springer Science & Business Media, 2009. [[pdf](https://sites.math.washington.edu/~rtr/papers/rtr169-VarAnalysis-RockWets.pdf)]
* \biblabel{Correa92}{Correa et al. (1992)} Correa, Rafael, Alejandro Jofre, and Lionel Thibault. "Characterization of lower semicontinuous convex functions." Proceedings of the American Mathematical Society (1992): 67-72. [[pdf](https://www.jstor.org/stable/2159295?seq=1#metadata_info_tab_contents)]
* \biblabel{Poliquin00}{Poliquin et al. (2000)} Poliquin, René, R. Rockafellar, and Lionel Thibault. "Local differentiability of distance functions." *Transactions of the American mathematical Society* 352.11 (2000): 5231-5249. [[pdf](https://www.ams.org/journals/tran/2000-352-11/S0002-9947-00-02550-2/S0002-9947-00-02550-2.pdf)]

