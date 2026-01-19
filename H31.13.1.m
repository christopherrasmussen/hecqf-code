/* This example was provided by an anonymous referee. We are extremely grateful
   for their generous suggestions. Some comments here are original to the referee;
   others are added by the authors as they experiment. */

/* Code for H(31, 13).1, first listed curve, with a_i invariants
   [0, 0, 1, 186930*a - 427490, 58571989*a - 135261471]
   defined over K = QQ(a), with a^2 + a - 3 = 0. */

/* Construct K and the heavenly elliptic curve E over that field */

Delta_K := 13;
ell := 31;

/* Build K as in our paper. */

_<T> := PolynomialRing(Rationals());

if Delta_K mod 4 eq 1 then
    minpoly_K := T^2 - T - (Delta_K - 1)/4;
else
    minpoly_K := T^2 - (Delta_K/4);
end if;

K<a> := NumberField( minpoly_K );

coeffs := [0, 0, 1, 186930*a - 427490, 58571989*a - 135261471];
E := EllipticCurve(coeffs);

/* Compute a degree-(ell-1)/2 factor of the ell-th division polynomial */

KX<X> := PolynomialRing(K);
dp, _, _ := KX!DivisionPolynomial(E, ell);
fac := Factorisation(dp);

f1 := fac[1][1];

/* Unclear to me why this is guaranteed to be the right factor */

print "DivisionPolynomial Factor: ", f1;

/* Find a root of f1 over K(\zeta_{23}), which will give the x-coordinate of some
   point P inside E[23]. */

time Kell := Compositum(K, CyclotomicField(ell));
time Kell := OptimisedRepresentation(Kell);

time rts := Roots(f1, Kell);

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
P := ChangeRing(E, Kell)![x0, y0];

/* Print the point and then check it has order 23. */

print P;
print ell * P;

