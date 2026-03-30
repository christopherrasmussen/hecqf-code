/* MAGMA script to read in candidate curves and test if they are heavenly */

/* Block giving contact info and request for comments.

*/

// load in candidate curves as isogeny classes
load "heavenly-candidates-complete.m";

QT<T> := PolynomialRing(Rationals());

sorted_keys := Sort(SetToSequence(Keys(data)));
heavenly_count := 0;

for Dlnkey in sorted_keys do
    l := Dlnkey[1];
    l_string := IntegerToString(l);
    
    D := Dlnkey[2];
    D_string := IntegerToString(D);
    
    n := Dlnkey[3];
    n_string := IntegerToString(n);

    isog_label := "C(" cat l_string cat ", " cat D_string cat ")." cat n_string;
    printf ">>>>> Candidate class %o:    [%o curves total] \n", isog_label, #data[Dlnkey][2];

    /* Construct minpoly and then the number field K = Q(a)... */
    minpoly_coeffs := data[Dlnkey][1];
    minpoly_K := QT ! minpoly_coeffs;
    K<a> := NumberField(minpoly_K);

    /* MAGMA will express K elements in Kmu_l_rel_extn whether we 
       want this behavior or not. We use a dummy ring with variable 'b'
       just to get human-readable output. */
    Dummy<b> := PolynomialRing(Integers());

    /* Since D and l are both currently fixed, we can go ahead and compute the 
       composite extension K(mu_\ell) */
    
    Qmu_l<z> := CyclotomicField(l);
    Kmu_l<xi0> := Compositum(K, Qmu_l);

    /* We need Kmu_l as a relative extension over Qmu_l, 
       except when l = 2. This is because MAGMA seems to 
       do better checking whether an element is a square
       in Kmu_l when it is working with Kmu_l as a 
       relative field extension over Qmu_l. */

    if l gt 2 then
        Kmu_l_rel_extn<xi> := RelativeField( Qmu_l, Kmu_l );
    else
        Kmu_l_rel_extn<xi> := Kmu_l;
    end if;

    /* We want a polynomial ring over K available */
    KX<X> := PolynomialRing(K);

    printf "   >> Curves defined over K = Q(b), where b satisfies %o.\n", minpoly_K;
    printf "\n";

    ec_list := data[Dlnkey][2];

    for raw_ainvs in ec_list do
        ainvs := [];
        dummy_ainvs_string := "";
        for raw_aj in raw_ainvs do
            Append(~ainvs, K ! raw_aj);
            dummy_ainvs_string := dummy_ainvs_string cat Sprint(Dummy ! raw_aj) cat ", ";
        end for;
        E := EllipticCurve(ainvs);

        Delta_E := Discriminant(E);
        a_invs := aInvariants(E);

        prime_support := PrimeDivisors( Integers() ! Norm(Discriminant(E)));
        /* if #prime_support gt 1 then
            printf "*** Curve has bad reduction at %o \n", prime_support;
            print "    Investigate further.";
        end if; */

        printf "Checking curve %o with bad reduction at %o \n", dummy_ainvs_string, prime_support;

        Phi := DivisionPolynomial(E, l);
        Phi_all_factors := Factorisation(Phi);
        Phi_factor := Phi_all_factors[1][1];
        if Degree(Phi_factor) le 2 then
            leading_terms := Phi_factor;
        else
            leading_terms := Phi_factor - (Phi_factor mod X^(Degree(Phi_factor)-2));
        end if;

        /* Factorisation appears to sort factors by increasing degree
           but if there is a bottleneck it may be worth it to ensure
           Phi_factor is chosen of minimum possible degree */

        /* We look for roots of Phi_factor over K(mu_l) */
        div_poly_roots := Roots(Phi_factor, Kmu_l_rel_extn);
        
        if #div_poly_roots eq 0 then
            /* If Phi_factor has no roots, exit */
            print "!!!!!!!!!!!!!!!!!!!!!! No roots of Phi_factor found.";
        else
            /* we grab the x-coordinate of a possible torsion point */
            x0 := div_poly_roots[1][1];
            // print "Found a root x0 inside K(mu_l)";

            /* There is definitely an l-torsion point of the form
               (x0, y0) in E[l] with x0 in Kmu_l, we just don't
               know if y0 in Kmu_l. We grab the hyperelliptic
               polynomials pE and qE such that E has equation

               E: y^2 + qE(x) * y = pE(x) */
            
            pE, qE := HyperellipticPolynomials(E);

            /* Set q0 := qE(x0) and p0 := pE(x0) */

            q0 := Evaluate(qE, x0);
            p0 := Evaluate(pE, x0);

            /* Now y0 is a root of y^2 + q0*y - p0 == 0, and this
               quadratic polynomial has discriminant q0^2 + 4p0 */
            
            y0_disc := q0^2 + 4 * p0;

            /* we have a root y0 over K(mu_l) precisely if y0_disc is
               a square in K(mu_l)... */
            
            // y_disc_is_square, y_disc_sqrt := IsSquare( Kmu_l_rel_extn ! y0_disc );
            y_disc_is_square, y_disc_sqrt := IsSquare( y0_disc );

            if y_disc_is_square then
                /* It appears we have an l-torsion point over K(mu_l). We
                   build this point P0 = (x0, y0) and verify l * P0 == O on E. */

                y0 := (-q0 + y_disc_sqrt)/2;
                P0 := ChangeRing(E, Kmu_l_rel_extn) ! [x0, y0];
                lP0 := l * P0;

                if IsZero(lP0) then
                    print " ";
                    printf "    >> %o-torsion point found: P = %o", l, P0;
                    print " ";
                    printf "    >> E[%o](K(mu_%o)) is nontrivial hence E heavenly at l = %o.", l, l, l;
                    heavenly_count := heavenly_count + 1;
                else
                    printf "    >> Something is seriously wrong; %o*P0 != O...", l;
                end if;

            else
                /* y_disc was not a square, so no torsion point found. */
                printf "    >> E(K(mu_%o)) = {O}, E is not heavenly.", l;
            end if;
            print " ";
        end if;
    end for;
end for;

printf "Found %o distinct heavenly pairs (E, l) with E heavenly at l. \n", heavenly_count;
