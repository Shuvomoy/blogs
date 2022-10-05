@def title = "References for my poster at the Cornell Young Researchers Workshop 2022"
@def published ="Oct 3, 2022"
@def tags =["optimization"]


# Cornell Young Researchers Workshop 2022
**Shuvomoy Das Gupta**

*October 6-8, 2022*

You have landed on this page, because you have scanned the QR code of my poster for Cornell Young Researchers Workshop 2022 😃.

Hello! I am very grateful that my research interests you 😊! In this page, I want to briefly provide you with some resources regarding my poster. Please feel free to send me an email 📧 to [sdgupta@mit.edu](mailto:sdgupta@mit.edu) if you have any questions regarding the poster, or if you just want to say hi!

---

**Table of contents**

\toc

---

## BnB-PEP paper

The poster (pdf  [here](https://shuvomoy.github.io/assets/Sozi_presentations/BnB_PEP_poster_final.pdf)) is based on the following paper.

> **Shuvomoy Das Gupta**, Bart P.G. Van Parys, and Ernest K. Ryu, “[Branch-and-Bound Performance Estimation Programming: A Unified Methodology for Constructing Optimal Optimization Methods](https://arxiv.org/abs/2203.07305)”, 2022. [[pdf](https://optimization-online.org/wp-content/uploads/2022/03/8819.pdf)] [[code](/assets/Sozi_presentations/BnB_PEP_poster_final.pdf)]

A detailed YouTube video describing the paper is below.

~~~
<iframe width="500" height="300" src="https://www.youtube.com/embed/sdYYFRxqbKQ" title="BnB-PEP: A Unified Methodology for Constructing Optimal Optimization Methods" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
~~~

## Materials to learn about PEP

BnB-PEP builds on PEP (performance estimation problem). I am very excited about PEP in general, so I want to tell you more about it!

### Papers

The performance estimation methodology, initiated by Drori and Teboulle in [1], formulates the worst case-performance of an optimization method as an optimization problem itself and upper bounds this performance through a semidefinite program (SDP) "relaxation. [Taylor](https://adrientaylor.github.io/), Hendrickx, and Glineur then showed in [2] that the SDP formulation is tight (not a relaxation) through the notion of "convex interpolation". Taylor, Hendrickx, and Glineur subsequently extended PEP to composite convex optimization in [3]. There are many more fantastic papers on PEP, but we can start learning about PEP by studying those three first:

[1] Y. Drori and M. Teboulle. Performance of first-order methods for smooth convex minimization: A novel approach. Mathematical Programming, 145(1-2):451–482, 2014. (link: [https://arxiv.org/pdf/1206.3209.pdf] (https://arxiv.org/pdf/1206.3209.pdf)

[2] A. B. Taylor, J. M. Hendrickx, and F. Glineur. Smooth strongly convex interpolation and exact worst-case performance of first-order methods. Mathematical Programming, 161(1-2):307–345, 2017. (link: [https://arxiv.org/pdf/1502.05666.pdf](https://arxiv.org/pdf/1502.05666.pdf)) 

[3] A. B. Taylor, J. M. Hendrickx, and F. Glineur. Exact worst-case performance of first-order methods for composite convex optimization. SIAM Journal on Optimization, 27(3):1283–1313, 2017. (link: [https://arxiv.org/pdf/1512.07516.pdf](https://arxiv.org/pdf/1512.07516.pdf))

### Blog

If you are interested to learn about PEP more informally without reading the papers above first, a gentle introduction is the following blog post by [Adrien Taylor](https://adrientaylor.github.io/):

[1] [A. B. Taylor, Computer-aided analyses in optimization](https://francisbach.com/computer-aided-analyses/)

### Video presentation

After the blog, a very nice video presentation by Adrien Taylor is available at this [link](https://stream.univie.ac.at/media/mathematik/owos/Adrien_Taylor-OWOS_Talk?res=1936). The slides for his talk is available [here](https://owos.univie.ac.at/fileadmin/user_upload/k_owos/Adrien_Taylor-OWOS.pdf).

A more advanced talk by Adrien Taylor on constructing optimal first-order optimization methods in convex optimization using PEP is below.

~~~
<iframe width="500" height="300" src="https://www.youtube.com/embed/jDgPIklp698" title="A. Taylor. A few constructive approaches to optimal first-order optimization methods in cvx. opt." frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
~~~







