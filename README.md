# hecqf-code

Code accompanying the paper [MR24]:

> Cam McLeman and Christopher Rasmussen, "Heavenly elliptic curves over quadratic fields,"
> [arXiv:2410.18389](https://arxiv.org/abs/2410.18389).

## Overview

This repository contains scripts for three calculations that support results appearing in [MR24].

* `balanced-filter.sage`: The balanced bound B(n,g) is shown to satisfy $B(2,1) < 1.68 \times 10^{7}$ in Corollary 5.3 of [MR24], but may be reduced to $B(2,1) \leq 19$ by a computational verification. The existence of a heavenly elliptic curve $E/K$ which is not balanced implies the Tate-Oort numbers $(j_{1}, j_{2})$ satisfy $j_{1} + j_{2} = e$ and the congruence (5.5):
$$M := \tau_{e} - q^{j_{1}} - q^{j_{2}} \equiv 0 \pmod{\ell}$$
Here, $\mathfrak{p}$ is a prime of $K$ dividing $p$, a rational prime, $q = \mathbf{N}\mathfrak{p} = p^{f}$, and $\tau_{e}$ is the trace of Frobenius for $\mathfrak{p}^{e}$. The crucial observation is that (5.5) must hold for every $\mathfrak{p} \nmid \ell$. 
Taking $p = 2$ and factoring $M$ for every possible choice of $f$ and $\tau_{e}$, we obtain a finite set of possible *unbalanced* $\ell$. The script then checks the congruence (5.5) against primes $2 < p \leq 11$, removing any primes $\ell$ that fail to satisfy (5.5) for all $p$ *with the same $(j_{1}, j_{2})$ for each $p$*.

* `cm-heavenly-search.sage`: A SageMath script that exhaustively searches for candidate heavenly elliptic curves. The output is a provably complete list of triples $(K, [E]_{K}, \ell)$ that necessarily contains all possible heavenly CM curves (up to $K$-isomorphism). As explained in [MR24], curves with $\mathbf{Q}$-rational $j$-invariant are excluded. A detailed explanation of the algorithm is given in [MR24, Sec. 7], as well as in the comments of the `.sage` file itself. Curves found with this program are only *candidates*. They may or may not actually be heavenly; this question is answered by the third script. Candidate curves are recorded in files `heavenly-candidates-complete.sage` and `heavenly-candidates-complete.m` for further use in `SageMath` or `MAGMA`, respectively.
  
* `heavenly-test.m`: A `MAGMA` script that reads in an array of isogeny classes of possibly-heavenly CM elliptic curves, and then tests whether the curve is heavenly by checking the existence or non-existence of an $\ell$-torsion point over the field $K(\mu_{\ell})$.  It suffices to check this for only one curve in each isogeny class to get a provably complete result; however, the calculation is sufficiently fast that we check every curve in every isogeny class for redundancy. (We thank an anonymous referee for outlining how to approach this calculation in `MAGMA`.)

## Usage

The balanced filter runs in a few seconds:

```bash
sage balanced-filter.sage
```

The CM heavenly search is computationally intensive (approximately 1.5 -- 2 hours on a basic paid `CoCalc` account as of Spring 2026); this is largely due to the sorting of curves into isogeny classes. It is recommended to redirect output to a log file; see the comments in the file for details.

```bash
sage cm-heavenly-search.sage
```

The Magma verification script should be run after the search completes:

```bash
magma -b heavenly-test.m
```

By default, `heavenly-test.m` merely reports the existence of an appropriate $\ell$-torsion point; it does not print the coordinates of such a point to `stdout`. To have the script print the points explicitly, use:

```bash
magma -b show_torsion:=y heavenly-test.m
```


## Authors

* Cam McLeman (University of Michigan--Flint)
* Christopher Rasmussen (Wesleyan University)

Comments welcome: crasmussen at wesleyan dot edu

## References

[BCP97] W. Bosma, J. Cannon, and C. Playoust, "The Magma algebra system. I.
The user language," *Journal of Symbolic Computation* **24** (1997),
pp. 235--265. Available at [http://magma.maths.usyd.edu.au/](http://magma.maths.usyd.edu.au/)

[Deu57] M. Deuring, "Die Zetafunktion einer algebraischen Kurve vom Geschlechte
Eins (iv)," *Nachrichten der Akademie der Wissenschaften in Göttingen*,
1957, pp. 55--80.

[DL15] H. Daniels and Á. Lozano-Robledo, "On the number of isomorphism classes of CM elliptic curves defined over a number field," *J. Number Theory* **157** (2015), pp. 367--396. 

[K02] M. Kida, "Potential good reduction of elliptic curves," *J. Symbolic Computation* (2002) **34**, pp. 173--180. [doi:10.1006/jsco.2002.0555](https://doi.org/10.1006/jsco.2002.0555)

[Lan73] S. Lang, *Elliptic Functions*, Addison-Wesley, 1973.

[MR24] C. McLeman and C. Rasmussen, "Heavenly elliptic curves over quadratic
fields," preprint, 2024. [arXiv:2410.18389](https://arxiv.org/abs/2410.18389).

[Sage] The Sage Developers, *SageMath, the Sage Mathematics Software System*.
Available at [https://www.sagemath.org/](https://www.sagemath.org/)
