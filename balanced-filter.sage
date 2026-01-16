r"""
Code to screen for cases where a heavenly elliptic curve could be *not* balanced, specifically in the
case where E is defined over a quadratic number field.

REFERENCE: "Heavenly elliptic curves over quadratic fields," 
           https://arxiv.org/pdf/2410.18389

AUTHORS: Cam McLeman (University of Michigan-Flint) and Christopher Rasmussen (Wesleyan University)

Comments welcome: crasmussen 'typical email symbol' wesleyan 'typical email punctuation' edu
"""

# We hard-code the construction of possible Tate-Oort numbers (aka "j-vectors") for practical use:

poss_jvec = []
for e in [1, 2, 3, 4, 6, 8, 12]:
    j1 = 0
    while (j1 < e/2):
        j2 = e - j1
        poss_jvec.append( (j1, j2) )
        j1 += 1

def possible_frob_trace(p, f):
    r'''Returns a list of all possible values for the trace of a Frobenius
        element for a prime pp of norm p^f. This is precisely

        ZZ \cap [-2\sqrt{p^f}, 2\sqrt{p^f}].'''

    q = p ** f
    extreme_frob_trace = 2 * q.isqrt()
    return list(range(-extreme_frob_trace, extreme_frob_trace + 1))

def trace_frob_power( trace_frob, pow, prime_norm ):
    r'''Suppose theta is a Frobenius element for a prime pp of
        norm prime_norm, and that P(T) is the integer characteristic polynomial
        for theta. Let alpha, beta be the complex roots of P(T), so that

        trace_frob = alpha + beta.

        This function returns the integral value of the trace of theta^pow,
        which is necessarily

        alpha^pow + beta^pow.

        The approach is to work recursively with symmetric polynomials.'''

    # we use two recursions. First, if pow is even, we can work with pow/2.
    if pow % 2 == 0:
        # We use the identity (a^{2m} + b^{2m}) = (a^m + b^m)^2 - 2(ab)^m
        half_pow = pow/2
        return trace_frob_power( trace_frob, half_pow, prime_norm )**2 - 2 * prime_norm ** half_pow
    # otherwise, use a simpler recursion
    elif pow > 1:
        # We use the identity a^m + b^m = (a^{m-1} + b^{m-1})(a + b) - ab(a^{m-2} + b^{m-2})
        return trace_frob_power( trace_frob, pow - 1, prime_norm ) * trace_frob - prime_norm * trace_frob_power( trace_frob, pow - 2, prime_norm)
    elif pow == 1:
        return trace_frob
    else:
        raise Exception("pow must be a positive integer.")

def alt_trace_frob_power( trace_frob, pow, prime_norm ):
    r'''An alternative to the function trace_frob_power(), which instead calculates the roots of the
        characteristic polynomial explicitly and then calculates the trace of the pow-th power directly.

        NOTE: When the characteristic polynomial is irreducible, this function is much slower than
        trace_frob_power(), because a number field construction is required.'''
    QQT.<T> = QQ[]
    char_poly = T^2 - ZZ(trace_frob) * T + ZZ(prime_norm)
    if char_poly.discriminant() == 0:
        # char poly splits over QQ, repeated root
        alpha = char_poly.roots()[0][0]
        return alpha ** pow + alpha ** pow
    elif char_poly.discriminant().is_square():
        # char poly factors over QQ with distinct roots
        alpha = char_poly.roots()[0][0]
        beta = char_poly.roots()[1][0]
        return alpha ** pow + beta ** pow
    else:
        # char poly is irreducible.
        K.<xi> = NumberField(char_poly)
        KT.<T> = K[]
        char_poly_over_K = KT(char_poly)
        alpha = char_poly_over_K.roots()[0][0]
        beta = char_poly_over_K.roots()[1][0]
        return ZZ(alpha ** pow + beta ** pow)
    raise ValueError("You shouldn't be here; check input.")

def tate_oort_trace( p, f, jvec ):
    r'''Calculates the expected trace modulo ell for the e-th power of Frobenius when
        the Tate-Oort numbers are in jvec. For a prime of norm q = p^f, this is the expression q^j1 + q^j2,
        which (for a heavenly elliptic curve) must be congruent to the trace of the e-th power of Frobenius, modulo ell.'''
    q = p ** f
    if len(jvec) != 2:
        raise Exception("jvec is not an acceptable vector of Tate-Oort numbers.")
    j1, j2 = jvec
    if j1 not in ZZ or j2 not in ZZ or j1 < 0 or j2 < 0:
        raise Exception("Tate-Oort numbers must be nonnegative integers.")
    return (q ** j1 + q ** j2)

def initial_screen( jvec, p0 ):
    r'''Given Tate-Oort numbers (j1, j2), calculates a list of primes ell satisfying:
           * ell > 11
           * gcd(e, ell - 1) divides j1   (e = j1 + j2)
           * there exists a prime pp of norm p0 or p0**2 such that
               tau_e - q**j1 - q**j2 % ell == 0,
             where q is the prime norm (i.e., p0 or p0**2), and tau_e is a possible integer
             trace value for the eth power of frobenius at pp.
    '''

    poss_ell = []
    e = sum( jvec )
    for inertial_degree in [1, 2]:
        prime_norm = p0 ** inertial_degree
        for tau in possible_frob_trace( p0, inertial_degree ):
            # The following if block is solely meant as a safety check for the trace of frobenius powers; one
            # could comment out this block and likely get a much faster execution time. Since
            # execution time is not an issue in this use case (abelian variety of dimension 1, field of definition of degree 2),
            # we leave the safety check on.
            if trace_frob_power( tau, e, prime_norm) != alt_trace_frob_power( tau, e, prime_norm ):
                print("ERROR CASE:", tau, e, prime_norm)
                print("trace: ", trace_frob_power( tau, e, prime_norm), "\t alt_trace: ", alt_trace_frob_power( tau, e, prime_norm ))
                print("jvec, p0, f: ", jvec, p0, inertial_degree)
            tau_e = trace_frob_power( tau, e, prime_norm )
            special_qty = tau_e - prime_norm ** jvec[0] - prime_norm ** jvec[1]
            # The special_qty must be zero mod ell, so we collect its prime divisors,
            # but we only keep ell > 11 that satisfy gcd(e, ell - 1) | j1.
            poss_ell_for_tau = [ell for ell in prime_divisors(special_qty) if (ell > 11 and jvec[0] % gcd(e, ell - 1) == 0)]
            for ell in poss_ell_for_tau:
                if ell not in poss_ell:
                    poss_ell.append(ell)
    poss_ell.sort()
    return poss_ell

# Using the functions above, we build a dictionary whose keys are Tate-Oort numbers (j1, j2)
# and whose values are lists of possible primes ell which might admit a heavenly elliptic curve
# with unbalanced Tate-Oort numbers (j1, j2)

poss_ell_dict = {}
for jvec in poss_jvec:
    list_of_ell_for_jvec = initial_screen( jvec, 2 )
    if len(list_of_ell_for_jvec) > 0:
        poss_ell_dict[ jvec ] = list_of_ell_for_jvec

# The next step will be to eliminate most of the remaining possibilites by checking constraints involving
# other small primes p0. Since p0 != ell is required, we (continue to) assume ell > 11 and use only
# the primes p0 with 2 < p0 <= 11. 

def congruence_check( jvec, ell, p0 ):
    r'''Given a vector jvec of Tate-Oort numbers, a prime ell, and an auxiliary prime p0 != ell,
        let theta be a Frobenius for a prime above p0. This function checks whether there is any
        possible choice for tau, the trace of theta, and f, the inertia degree of a prime over p0,
        such that the congruence 

            tau_e - q ** j1 - q ** j2 == 0 (modulo ell) 

        holds. Here, e = j1 + j2, and tau_e is the trace of theta^e.

        Returns True if some choice exists; returns False if no choice exists.'''

    if p0 == ell:
        raise ValueError("The prime p0 must be distinct from ell.")

    e = sum( jvec )
    for inertial_degree in [1, 2]:
        prime_norm = p0 ** inertial_degree
        for tau in possible_frob_trace( p0, inertial_degree ):
            tau_e = trace_frob_power( tau, e, prime_norm )
            special_qty = tau_e - prime_norm ** jvec[0] - prime_norm ** jvec[1]
            if special_qty % ell == 0:
                return True
    return False

for p0 in [3, 5, 7, 11]:
    remaining_jvec = list(poss_ell_dict.keys())
    for jvec in remaining_jvec:
        new_list = []
        for ell in poss_ell_dict[ jvec ]:
            if congruence_check( jvec, ell, p0 ):
                new_list.append( ell )
        if len(new_list) > 0:
            poss_ell_dict[ jvec ] = new_list
        else:
            del poss_ell_dict[ jvec ]

# Note that poss_ell_dict drops a key jvec if there are no possible ell > 11 left for jvec.
#
# Thus, the final status of poss_ell_dict is that its keys are Tate-Oort numbers that might
# admit an unbalanced ell > 11, and its value is a list of those possible ell.
#
# Report the results to the user:

print("The remaining unbalanced cases with ell > 11 are as follows. \n")

print("Tate-Oort Numbers  Possible ell > 11")
print("-----------------  -----------------")   
for jvec in poss_ell_dict:
    print( str(jvec).rjust(13), "      ", poss_ell_dict[ jvec ] )
