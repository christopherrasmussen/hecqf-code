/* Loads the heavenly candidates curves database curves_db from HeavenlyCandidates.txt */

load "curve_data.m";

CheckHeavenlyCurve := function(Delta_K, ell, raw_coeffs)

    /*  Expects a triple of data representing an elliptic curve heavenly over a quadratic field K 
         Delta_K represents the discriminant of K and ell the prime of heavenliness 
         raw_coeffs expects [a1, a2, a3, a4, a6] where each coefficient ai is a tuple [c0,c1]
         representing the point c0+c1*a where a is a specific field generator defined below 
    */

    print "--------------------------------------------------------";
    printf "Checking Curve: Delta_K = %o, ell = %o \n", Delta_K, ell;
    print "--------------------------------------------------------";

    /* 1. Construct K and E */
    _<T> := PolynomialRing(Rationals());

    if Delta_K mod 4 eq 1 then
        minpoly_K := T^2 - T - (Delta_K - 1)/4;
    else
        minpoly_K := T^2 - (Delta_K/4);
    end if;

    K<a> := NumberField( minpoly_K );
    // printf "Field K defined by: %o\n", minpoly_K;

    coeffs := [];
    for item in raw_coeffs do
        Append(~coeffs, K!(item[1] + item[2]*a));
    end for;

    E := EllipticCurve(coeffs);

    /* 2. Factor the Division Polynomial  */

    Qmu<xi> := CyclotomicField( ell );
    Kmu := Compositum(K, Qmu);
    
    KX<X> := PolynomialRing(K);
    dp := DivisionPolynomial(E, ell);
    fac := Factorisation(dp);
    
    f1 := fac[1][1];

    printf "Degree of Div Poly Factor: %o \n", Degree(f1);

    /* 3. Find x-coordinate for P */

    rts := Roots(f1, Kmu);
    if #rts eq 0 then
        print ">> FAILURE: No roots found for this factor in Kmu.";
        return false;
    end if;

    x0 := rts[1][1];
    print "Found x-coordinate of P.";

    /* 4. Find coordinate for P */
    qE, pE := HyperellipticPolynomials(E);
    // printf "coeffs: %o \n", coeffs;
    // printf "qE: %o \n", qE;
    RHS := Evaluate(qE, x0);
    LHS := Evaluate(pE, x0);
    y_disc := LHS^2 + 4 * RHS;

    // printf "Field Qmu is %o \n", Qmu;
    // printf "Field Kmu is %o \n", Kmu;

    /* Build relative field for faster square-checking */
    
    if ell gt 2 then
        Kmu_copy := RelativeField( Qmu, Kmu );
    else
        Kmu_copy := Kmu;
    end if;

    // printf "Field Kmu_copy is %o \n", Kmu_copy;

    is_sq, y_disc_root := IsSquare( Kmu_copy!y_disc );

    if not is_sq then
        second_try_is_sq, second_disc_root := IsSquare( y_disc );
        if not second_try_is_sq then
            print ">> FAILURE: y_disc is not a square in Kmu.";
            return false;
        end if;
    end if;

    y0 := (-LHS + y_disc_root)/2;
    coerced_x0 := Kmu_copy!x0; 

    /* 6. Verify Computation */
    P := ChangeRing(E, Kmu_copy)![coerced_x0, y0];
    print "Point P constructed.";
    check_point := ell * P;
    
    print "Verifying Order...";
    // print "P =", P;
    printf "ell * P = %o\n", check_point;
    
    if IsZero(check_point) then
        print ">> SUCCESS: Point has correct order.";
        return true;
    else
        print ">> FAILURE: Point does NOT have correct order.";
        return false;
    end if;

end function;

for i in [1..#curves_db] do
    ell := curves_db[i][2];
    CheckHeavenlyCurve(curves_db[i][1], curves_db[i][2], curves_db[i][3]);
end for;