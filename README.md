# hecqf-code

Code accompanying the paper:

> Cam McLeman and Christopher Rasmussen, "Heavenly elliptic curves over quadratic fields,"
> [arXiv:2410.18389](https://arxiv.org/abs/2410.18389).

## Overview

This repository contains three scripts that together perform an exhaustive
search for CM elliptic curves over quadratic fields that are heavenly, and
verify the results. The computation proceeds in three stages:

1. **Balanced filter** (`balanced-filter.sage`): A SageMath script that screens
   for primes ell > 11 where an unbalanced heavenly curve could exist. It
   uses Frobenius trace congruences at small auxiliary primes to eliminate
   most possibilities.

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
