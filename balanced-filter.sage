r"""
Code to screen for cases where a heavenly elliptic curve could be *not* balanced, specifically in the
case where E is defined over a quadratic number field.

REFERENCE: "Heavenly elliptic curves over quadratic fields,"
           https://arxiv.org/pdf/2410.18389

AUTHORS: Cam McLeman (University of Michigan-Flint) and Christopher Rasmussen (Wesleyan University)

Comments welcome: crasmussen 'typical email symbol' wesleyan 'typical email punctuation' edu
"""

# We construct all possible "unbalanced" Tate-Oort numbers
# (j1, j2) satisfying j1 + j2 == e and 0 <= j1 < j2 < e/2.
# These are stored in poss_jvec.

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
    extreme_frob_trace = floor(2 * q.sqrt())
    return list(range(-extreme_frob_trace, extreme_frob_trace + 1))

def trace_frob_power( trace_frob, pow, prime_norm ):
    r'''Suppose theta is a Frobenius element for a prime pp of
        norm prime_norm, and that P(T) is the integer characteristic
        polynomial of theta. Let alpha, beta be the complex roots
        of P(T), so that

        trace_frob = alpha + beta.

        This function returns the integral value of the trace of theta^pow,
        which is necessarily

        alpha^pow + beta^pow.

        The approach is to work recursively with symmetric polynomials.'''

    # we use two recursions. First, if pow is even, we can work with pow/2.
    if pow % 2 == 0:
        # We use the identity (a^{2m} + b^{2m}) = (a^m + b^m)^2 - 2(ab)^m
        half_pow = pow // 2
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

        NOTE: When the characteristic polynomial is irreducible, this function is slower than
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
        KT.<T0> = K[]
        char_poly_over_K = KT(char_poly)
        alpha = char_poly_over_K.roots()[0][0]
        beta = char_poly_over_K.roots()[1][0]
        return ZZ(alpha ** pow + beta ** pow)

def initial_screen( jvec, p0 ):
    r'''Given Tate-Oort numbers (j1, j2), calculates a list of primes ell satisfying:
           * ell > 11
           * gcd(e, ell - 1) divides j1   (e = j1 + j2)
           * (tau_e - q ** j1 - q ** j2) % ell == 0 for *at least one*
             choice of q in [p, p**2] and some value tau_e
             which is a valid Frobenius trace for a (hypothetical)
             prime pp of norm q.
    '''

    poss_ell = []
    e = sum( jvec )
    for inertial_degree in [1, 2]:
        prime_norm = p0 ** inertial_degree
        for tau in possible_frob_trace( p0, inertial_degree ):
            # The following if block is solely meant as a safety check for the
            # trace of frobenius powers; one could comment out this block and
            # likely get a much faster execution time. Since execution time
            # is not an issue in this use case -- abelian variety of dimension 1,
            # field of definition of degree 2 -- we leave the safety check in place.
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

# The next step will be to eliminate most of the remaining possibilities by checking constraints involving
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

def compat_check_fixed_f_fixed_tau(ell, p0, f, ivec, jvec, tau):
    r'''Checks compatibility between a pair of heavenly exponents ivec and
        a pair of Tate-Oort numbers jvec. Returns True if and only if:

          - e * i[r] == j[r] (mod (ell-1)) for r = 0, 1,
          - tau_e == q0 ** j[0] + q0 ** j[1] (mod ell),
          -   tau == q0 ** i[0] + q0 ** i[1] (mod ell).'''
    q0 = p0 ** f
    e = sum(jvec)
    for r in [0, 1]:
        if (e*ivec[r] - jvec[r]) % (ell-1) != 0:
            return False
    tau_e = alt_trace_frob_power(tau, e, q0)
    if (tau_e - q0 ** jvec[0] - q0 ** jvec[1]) % ell != 0:
            return False
    if (tau - q0 ** ivec[0] - q0 ** ivec[1]) % ell != 0:
            return False
    return True

def compat_check_fixed_f(ell, p0, f, ivec, jvec):
    r'''Returns True if there exists at least one tau for which
        compat_check_fixed_f_fixed_tau returns True. Otherwise returns False.'''
    for tau in possible_frob_trace(p0, f):
        if compat_check_fixed_f_fixed_tau(ell, p0, f, ivec, jvec, tau):
            return True
    return False

def compat_check(ell, p0, ivec, jvec):
    r'''Returns True if there exists f in [1, 2] for which
        compat_check_fixed_f returns True. Otherwise returnse False.'''
    for f in [1, 2]:
        if compat_check_fixed_f(ell, p0, f, ivec, jvec):
            return True
    return False
            
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

# We are not done. Suppose pp0 is a prime ideal of K lying above the rational prime p0
# with norm q0 = p0 ** f and p0 != ell. Let theta0 be a Frobenius element for pp0. Because
# the hypothetical elliptic curve E/K is heavenly, it has an l-torsion Galois representation
# of the form:
#
#    rho_{E,ell} : G_K --> [ \chi^{i1}     *     ]
#                          [     0     \chi^{i2} ]
#
# where \chi is the mod-l cyclotomic character. For E/K to exist there must be a fixed choice
# of Tate-Oort numbers (j1, j2) and a fixed pair of exponents (i1, i2) that satisfy all of
# the following conditions:
#
#  (1)    e * ir == jr (mod (ell - 1))    for r = 1, 2
#
#  (2)    For every p0 != ell, there exists an f in [1, 2] and an integer tau (a possible Frobenius trace at p0) such that:
#
#         (a) |tau| \leq 2 * sqrt(p0 ** f)
#         (b)    tau == (p0^f)^i1 + (p0^f)^i2   (mod ell)
#         (c)  tau_e == (p0^f)^j1 + (p0^f)^j2   (mod ell)
#
#         Here, tau_e is the trace of the eth power of the same Frobenius element.
#
#  The point is that if (2) fails, then it is impossible for rho_{E,ell} to "behave correctly" if E is heavenly with
#  Tate-Oort numbers j1, j2.
#
#  So we screen the remaining cases of (ell, (j1, j2)) against possible (i1, i2) by brute force.
#
#  A table is created for each remaining ell. Each row is indexed by a not-balanced (j1, j2)
#  Each column is indexed by i1 for a possible pair of exponents (i1, i2). (Note that fixing i1 determines i2.)
#  The entry in the cell for (j1, j2) and (i1, i2) is:
#
#  --      if the congruence (1) fails
#  p0      if (2) fails for prime p0
#  blank   if this pair (j1, j2) remains viable for a heavenly E/K.

print(" ")
print("Screening remaining unbalanced cases against individual (i1, i2)... ")
print(" ")

remaining_unbalanced_cases = {}

remaining_cases = []
for jvec in poss_ell_dict:
    for ell in poss_ell_dict[jvec]:
        if ell in remaining_unbalanced_cases:
            remaining_unbalanced_cases[ ell ].append( jvec )
        else:
            remaining_unbalanced_cases[ ell ] = [ jvec ]

for ell in remaining_unbalanced_cases:
    print(f"*** ell = {ell} *** |   i1 of exponents (i1, i2) ")
    table_header = "* ( j1, j2) *    | "
    for i1 in [0..(ell-2)]:
        table_header += " " + str(i1).rjust(2) + " "
    print(table_header)
    separator = (11 + 4*(ell-1) + 7)*"-"
    print(separator)
    for jvec in remaining_unbalanced_cases[ ell ]:
        eliminated_flag = True
        e = sum(jvec)
        j1, j2 = jvec
        table_row = "  ( "+str(j1).rjust(2)+", "+str(j2).rjust(2)+")      | "
        for i1 in [0..ell-2]:
            i2 = ell - i1
            ivec = (i1, i2)
            if (e*i1 - j1) % (ell-1) != 0:
                table_row += " -- "
            elif not compat_check(ell, 2, ivec, jvec):
                table_row += "  2 "
            elif not compat_check(ell, 3, ivec, jvec):
                table_row += "  3 "
            elif not compat_check(ell, 5, ivec, jvec):
                table_row += "  5 "
            else:
                table_row += "    "
                eliminated_flag = False
        if eliminated_flag:
            table_row += "Impossible"
        print(table_row)
    print(" ")
