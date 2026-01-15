"""
Code to screen for cases where a heavenly elliptic curve could be unbalanced, specifically in the
case where E is defined over a quadratic number field.
"""

# package imports

# we hard-code the construction of possible tate-oort numbers (aka "j-vectors") for practical use
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
        norm prime_norm, and that P(T) is the integer characteristic polyomial
        for theta. Let alpha, beta be the complex roots of P(T), so that

        trace_frob = alpha + beta.

        This function returns the integral value of the trace of theta^pow,
        which is necessarily
        
        alpha^pow + beta^pow.'''
    
    # we use two recursions. First, if pow is even, we can take work with pow/2.
    if pow % 2 == 0:
        half_pow = pow/2
        return trace_frob_power( trace_frob, half_pow, prime_norm )**2 - 2 * prime_norm ** half_pow
    # otherwise, use a simpler recursion
    elif pow > 1:
        return trace_frob_power( trace_frob, pow - 1, prime_norm ) * trace_frob - prime_norm * trace_frob_power( trace_frob, pow - 2, prime_norm)
    elif pow == 1:
        return trace_frob
    else:
        raise Exception("pow must be a positive integer.")

def tate_oort_trace( p, f, jvec ):
    r'''Calculates the expected trace modulo ell for the e-th power of Frobenius when
        the Tate-Oort numbers are in jvec. For a prime of norm q = p^f, this is the expression q^j1 + q^j2,
        which we require to be congruent to the trace of the e-th power modulo ell.'''
    q = p ** f
    if len(jvec) != 2:
        raise Exception("jvec is not an acceptable vector of Tate-Oort numbers.")
    j1, j2 = jvec
    if j1 not in ZZ or j2 not in ZZ or j1 < 0 or j2 < 0:
        raise Exception("Tate-Oort numbers must be nonnegative integers.")
    return (q ** j1 + q** j2)

def initial_screen( jvec, p0 ):
    r'''Given Tate-Oort numbers (j1, j2), calculates a list of primes ell satisfying:
           * ell > 11
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
            tau_e = trace_frob_power( tau, e, prime_norm )
            special_qty = tau_e - prime_norm ** jvec[0] - prime_norm ** jvec[1]
            poss_ell_for_tau = [ell for ell in prime_divisors(special_qty) if ell > 11]
            for ell in poss_ell_for_tau:
                if ell not in poss_ell:
                    poss_ell.append(ell)
    poss_ell.sort()
    return poss_ell

# using the functions above, we build a dictionary whose keys are Tate-Oort numbers (j1, j2)
# and whose values are lists of possible primes ell which might admit a heavenly elliptic EllipticCurve
# with unbalanced Tate-Oort numbers (j1, j2)

poss_ell_dict = {}
for jvec in poss_jvec:
    poss_ell_dict[ jvec ] = initial_screen( jvec, 3 )

