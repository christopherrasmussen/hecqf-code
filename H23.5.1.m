/* This example was provided by an anonymous referee. We are extremely grateful
   for their generous suggestions. Some comments here are original to the referee;
   others are added by the authors as they experiment. */

/* Code for H(23, 5).1, first listed curve, with a_i invariants
   [0, 0, 1, 46 * a - 368, 2645 * a - 6216]
   defined over K = QQ(a), with a^2 + a - 1 = 0. */

/* Construct K = Q(sqrt{5}) and the heavenly elliptic curve E over that field */

_<x> := PolynomialRing(Rationals());
K<a> := NumberField(x^2 + x - 1);
coeffs := [0, 0, 1, -46*a - 368, -2645*a - 6216];
E := EllipticCurve(coeffs);

/* Compute a degree-11 factor of the 23rd division polynomial */

KX<X> := PolynomialRing(K);
dp, _, _ := KX!DivisionPolynomial(E, 23);
fac := Factorisation(dp);


f1 := fac[1][1];

/* Unclear to me why this is guaranteed to be the right factor */

print "DivisionPolynomial Factor: ", f1;

/* Find a root of f1 over K(\zeta_{23}), which will give the x-coordinate of some
   point P inside E[23]. */

time K23 := Compositum(K, CyclotomicField(23));
time K23 := OptimisedRepresentation(K23);

time rts := Roots(f1, K23);

x0 := rts[1][1];   // the x-coordinate of P

/* Construct the corresponding y-coordinate */

hE := HyperellipticPolynomials(E);
LHS := Evaluate(hE, x0);

/* Check that the discriminant of y^2 + y = LHS is a square, and use the square root
   to construct the y-coordinate.

   Note that the equation of E is of the form y^2 + y = f(x). */

disc := 1 + 4 * LHS;
_, sqrDisc := IsSquare(disc);

y0 := (-1 + sqrDisc)/2;
P := ChangeRing(E, K23)![x0, y0];

/* Print the point and then check it has order 23. */

print P;
print 23*P;

