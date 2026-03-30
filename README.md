# hecqf-code

Code accompanying the paper [MR24]:

> Cam McLeman and Christopher Rasmussen, "Heavenly elliptic curves over quadratic fields,"
> [arXiv:2410.18389](https://arxiv.org/abs/2410.18389).

## Overview

This repository contains scripts for three calculations that support results appearing in [MR24].

1. `balanced-filter.sage`: The balanced bound B(n,g) is shown to satisfy $B(2,1) < 1.68 \times 10^{7}$ in Corollary 5.3 of [MR24], but may be reduced to $B(2,1) \leq 19$ by a computational verification. The existence of a heavenly elliptic curve $E/K$ which is not balanced implies the Tate-Oort numbers $(j_{1}, j_{2})$ satisfy $j_{1} + j_{2} = e$ and the congruence (5.5):
$$M := \tau_{e} - q^{j_{1}} - q^{j_{2}} \equiv 0 \pmod{\ell}$$
Here, $\mathfrak{p}$ is a prime of $K$ dividing $p$, a rational prime, $q = \mathbf{N}\mathfrak{p} = p^{f}$, and $\tau_{e}$ is the trace of Frobenius for $\mathfrak{p}^{e}$. The crucial observation is that (5.5) must hold for every $\mathfrak{p} \nmid \ell$. 
Taking $p = 2$ and factoring $M$ for every possible choice of $f$ and $\tau_{e}$, we obtain a finite set of possible *unbalanced* $\ell$. The script then checks the congruence (5.5) against primes $2 < p \leq 11$, removing any primes $\ell$ that fail to satisfy (5.5) for all $p$ *with the same $(j_{1}, j_{2})$ for each $p$*.

2. **CM heavenly search** (`cm-heavenly-search.sage`): A SageMath script that
   exhaustively searches for candidate heavenly elliptic curves. Starting
   from the finitely many CM j-invariants of class number 2, it identifies
   curves with good reduction away from a single prime (up to quadratic
   twist), then sorts them into isogeny classes. Output is written to
   `heavenly-candidates-complete.sage` and `heavenly-candidates-complete.m`.

3. **Heavenly test** (`heavenly-test.m`): A Magma script that reads the
   candidate curves produced by Step 2 and tests whether each admits
   nontrivial ell-torsion over the cyclotomic extension K(mu_ell), which is
   the defining condition for being heavenly.

## Usage

The balanced filter runs in a few seconds:

```
sage balanced-filter.sage
```

The CM heavenly search is computationally intensive (approximately 1.5--2
hours on a basic paid CoCalc account as of Spring 2026). It is recommended
to redirect output to a log file; see the comments in the file for details.

```
sage cm-heavenly-search.sage
```

The Magma verification script should be run after the search completes:

```
magma heavenly-test.m
```

## Authors

- Cam McLeman (University of Michigan--Flint)
- Christopher Rasmussen (Wesleyan University)

Comments welcome: crasmussen at wesleyan dot edu

## References

[BCP97] W. Bosma, J. Cannon, and C. Playoust, "The Magma algebra system. I.
The user language," *Journal of Symbolic Computation* **24** (1997),
pp. 235--265. Available at http://magma.maths.usyd.edu.au/

[Deu57] M. Deuring, "Die Zetafunktion einer algebraischen Kurve vom Geschlechte
Eins (iv)," *Nachrichten der Akademie der Wissenschaften in Göttingen*,
1957, pp. 55--80.

[Lan73] S. Lang, *Elliptic Functions*, Addison-Wesley, 1973.

[MR24] C. McLeman and C. Rasmussen, "Heavenly elliptic curves over quadratic
fields," preprint, 2024. [arXiv:2410.18389](https://arxiv.org/abs/2410.18389).

[Sage] The Sage Developers, *SageMath, the Sage Mathematics Software System*.
Available at https://www.sagemath.org/
