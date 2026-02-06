/* Code for H(47, 5).1, first listed curve, with a_i invariants
   [0, 0, 1, 4136 * a - 17578, 324723 * a - 962572]
   defined over K = QQ(a), with a^2 + a - 1 = 0. */

/* Construct K and the heavenly elliptic curve E over that field */

Delta_K := 5;
ell := 47;

/* Build K as in our paper. */

_<T> := PolynomialRing(Rationals());

if Delta_K mod 4 eq 1 then
    minpoly_K := T^2 - T - (Delta_K - 1)/4;
else
    minpoly_K := T^2 - (Delta_K/4);
end if;

K<a> := NumberField( minpoly_K );
Qmu<xi> := CyclotomicField( ell );
Kmu := Compositum(K, Qmu);

coeffs := [0, 0, 1, 4136*a - 17578, 324723 * a - 962572];
E := EllipticCurve(coeffs);

/* Compute a degree-(ell-1)/2 factor of the ell-th division polynomial */

KX<X> := PolynomialRing(K);
dp, _, _ := KX!DivisionPolynomial(E, ell);
fac := Factorisation(dp);

f1 := fac[1][1];

/* Unclear to me why this is guaranteed to be the right factor */

print "DivisionPolynomial Factor: ", f1;

time rts := Roots(f1, Kmu);

print("Found a root.");

x0 := rts[1][1];   // the x-coordinate of P

print("Got the x coord.");
print("x-coord of P is");
print(x0);

/* Construct the corresponding y-coordinate */

hE := HyperellipticPolynomials(E);
RHS := Evaluate(hE, x0);
y_disc := 1 + 4 * RHS;

print("calculation towards y coord done");

/* We build a copy of Kmu which is a relative extension over Qmu */

Kmu_copy := RelativeField( Qmu, Kmu );

/* Check if y poly discriminant is a square inside Kmu */

_, y_disc_root := IsSquare( Kmu_copy!y_disc );

y0 := (-1 + y_disc_root)/2;  // should be y-coord of a torsion point.
coerced_x0 := Kmu_copy!x0; 

P := ChangeRing(E, Kmu_copy)![coerced_x0, y0];

print("Calculated P");

/* Print the point and then check it has order 23. */

print P;
print ell * P;

