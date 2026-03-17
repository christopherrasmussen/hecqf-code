r"""
Code to exhaustively search for triples (K, [E]_K, l), where:

    - l is a rational prime,
    - K is a quadratic field,
    - E/K is a CM elliptic curve, heavenly at l,
    - j(E) not in QQ.

REFERENCE: "Heavenly elliptic curves over quadratic fields,"
           https://arxiv.org/pdf/2410.18389

AUTHORS: Cam McLeman (University of Michigan-Flint) and Christopher Rasmussen (Wesleyan University)

Comments welcome: crasmussen 'typical email symbol' wesleyan 'typical email punctuation' edu
"""

import time
true_start = time.time()

T = polygen(ZZ,'T')
def preferred_char_poly(Delta, indeterminate = T):
    r'''
    For a positive real integer Delta, returns a standardized form
    for the minimal polynomial of the unique quadratic field K of absolute
    discriminant Delta, in the specified indeterminate.'''

    if Delta % 4 == 0:
        return indeterminate^2 - (Delta/4)
    elif Delta % 4 == 1:
        return indeterminate^2 - indeterminate - (Delta - 1)/4
    else:
        raise Exception("Delta is not the discriminant of a quadratic field.")

def reduce_discriminant_by_qtwist(E0, l):
    r'''

    INPUT: E0, l, where
        -- E0, an elliptic curve
        -- l, a rational prime number

    OUTPUT: (twist_flag, E1), where
        -- twist_flag = True if a curve good outside {2, l} is found, or
           twist_flag = False if no such curve exists.
        -- E1, a quadratic twist of E0 with good reduction away
           from {2, l} (if twist_flag = True), or None (if twist_flag = False).'''

    E = E0.global_minimal_model()
    K = E.base_ring()
    OK = K.ring_of_integers()

    # We identify the primes we allow to appear in the minimal discriminant ideal.
    allowed_bad_primes = [ pp for pp, val in K.fractional_ideal(2*l).factor() ]

    # We initialize a constant which will eventually specify a quadratic twist.
    global_twist_constant = K(1)

    # We check each prime in the discriminant of E for a contribution to global_twist_constant

    # Note that the code relies on the fact that we will only see fields K with class number 1;
    # any generalization of this code may need to factor the minimal discriminant ideal instead.

    for pp, val in E.discriminant().factor():
        # We exclude primes pp that are allowed to be bad:
        if pp not in allowed_bad_primes:
            # We extract the largest sixth power possible.
            # To remove a^6, we twist by a^{-1}.
            local_twist_valuation = -(val // 6) # reminder, // is integral division w/o remainder
            global_twist_constant = global_twist_constant * (pp ** local_twist_valuation)

    E1 = E.quadratic_twist( global_twist_constant ).global_minimal_model()

    # now check if every remaining prime dividing the discriminant is allowed:
    for pp, val in E1.discriminant().factor():
        if pp not in allowed_bad_primes:
            twist_flag = False
            return (twist_flag, None)

    # every prime was allowed to be bad, so this candidate works
    twist_flag = True
    return (twist_flag, E1)

def quadratic_Sunit_list(E_input):
    # given an elliptic curve whose conductor has rational support over a set like {2, l},
    # where l is a prime divding 2m, and L = QQ(rt-m) is the CM field, find a list of
    # primes in K that we should consider twisting over
    cm_disc = E_input.cm_discriminant()
    K = E_input.base_ring()
    OK = K.ring_of_integers()
    NE_times_DeltaL = 2 * cm_disc * E_input.conductor()
    S_QQ = [p for p, val in NE_times_DeltaL.norm().factor()]
    S = []
    for p in S_QQ:
        for pp, val in (p*OK).factor():
            if pp not in S:
                S.append(pp)
    return S

def is_pair_isogenous(E1, l1, E2, l2):
    r'''
    INPUT: E1, l1, E2, L2, where
       -- E1 is an elliptic curve good outside l1
       -- E2 is an elliptic curve good outside l2
    OUTPUT: flag, where
       -- flag == True if the following three conditions are all met:
             -- l1 == l2
             -- E1.base_ring() is isomorphic to E2.base_ring()
             -- E1 is isogenous to E2 over the common base field
       -- flage == False in any other case.'''

    # check the primes of bad reduction
    if l1 != l2:
        return False
    l = l1

    # primes match, so compare base fields

    K1 = E1.base_ring()
    K2 = E2.base_ring()
    if not K1.is_isomorphic(K2):
        return False

    # fields of definition match agree. A necessary condition for E1/K and E2/K
    # to be K-isogenous is that the CM fields for E1 and E2 agree. This is a
    # much faster check than looking for an isogeny, so we screen against this
    # first.

    L1 = QuadraticField(E1.cm_discriminant())
    L2 = QuadraticField(E2.cm_discriminant())

    if not L1.is_isomorphic(L2):
        return False

    # Now that E1 and E2 have the CM field, all we can do is check for
    # an isogeny.

    if E1.is_isogenous(E2):
        return True
    else:
        return False

# The calculation, especially the sort into isogeny classes, is
# quite lengthy. Progress is written to the text file
#
# output.log
#
# in the present directory while running. If the calculation
# completes, results (all candidate CM curves that might
# be heavenly, up to isomorphism over their base fields)
# is written to the file
#
# heavenly-candidates-complete.sage
#
# also in the present directory. On a basic (paid)
# CoCalc account in Spring 2026, this calculation had
# a run-time of ~2 hours.
#
# Here is a brief overview of the strategy.
#
# 1. We consider each of the finitely many j-invariants which
#    correspond to an elliptic curve E/QQ(j(E)) with CM by an
#    order of class number 2. For each such j-invariant, we
#    construct an integral elliptic curve E0/QQ(j) with
#    j(E0) = j. We ignore QQ-rational j-invariants. This
#    generates a list of 58 curves.
#
# 2. We refer to Theorem 10 of Chapter 10, Section 4 of the
#    following reference:
#
#    Reference: Elliptic functions, Serge Lang,
#               Addison-Wesley, 1973.
#
#    though we note Lang attributes this result to Deuring:
#
#    Reference: Die Zetafunktion einer algebräischen Kurve
#               vom Geschlechte Eins, (iv), Nachrichten
#               Akad. Wiss. Göttingen, 1957, pp. 55--80.
#
#    We only need the first claim of Lang's statement of the
#    theorem. In the notation of our paper, the result can be
#    phrased as follows:
#
#    Let K be a number field, and suppose E/K is an elliptic
#    curve with complex multiplication by an order in the
#    imaginary quadratic field L. Suppose K is not contained
#    in L. If pp is a prime of K where E/K has good reduction,
#    then pp must be unramified in the extension KL/K.
#
#    (Because the field K = QQ(j(E)) must be real quadratic,
#    the hypothesis ``K not contained in L'' always holds.)
#
#    For each curve E/K from the list of 58 curves, we
#    calculate the ramification of KL/K, or rather its
#    support over QQ: That is, we specifically calculate
#    the following set:
#
#       RamQ(KL/K) = { p in ZZ : \exists pp in O_K with
#                      pp | p and pp ramified in KL/K }
#
#    We break into cases:
#
#    -- if # RamQ(KL/K) > 1, there is no possibility of
#       finding a heavenly curve that is a quadratic twist
#       of E/K; we abandon this curve.
#    -- if RamQ(KL/K) == {l}, then the *only* possible
#       prime at which a quadratic twist could be heavenly
#       is l. We proceed to Step 4 with the pair (E, l).
#    -- if RamQ(KL/K) == {}, then several heavenly primes
#       are still possible for some quadratic twist. We
#       proceed to Step 3.
#
# 3. Let Isog(E/K) be the set of rational primes l for which
#    E admits a K-rational isogeny of degree l. (From the
#    structure theorem for heavenly abelian varieties, we may
#    be sure that E has a K-rational isogeny of degree l if E
#    is heavenly at l.) It is known that Isog(E/K) is a finite
#    set. We calculate Isog(E/K) and then loop over the set:
#
#    -- For each l in Isog(E/K): Send (E, l) to Step 4.
#
# 4. Given (E,l), we calculate the minimal discriminant ideal
#    for E, Delta_{E/K}. For each prime ideal pp of O_K with
#    pp \nmid 2*l, we attempt to twist the factor of pp out
#    of Delta_{E/K}. Either this is successful for every such
#    pp, and we obtain E'/K, a twist over K of E/K with good
#    reduction outside {2,l}, or this fails for some pp, and
#    we abandon this pair. If we are successful, we pass
#    (E', l) to Step 5.
#
# 5. At this point we have an elliptic curve E'/K with good
#    reduction away from {2,l}, and must solve the following
#    problem: Find (up to K-isomorphism) all quadratic twists
#    E''/K of E'/K which have good reduction away from l. The
#    data of E''/K is completely determined by the quadratic
#    extension F/K which is used to twist E'/K. We use a
#    conductor relation for quadratic twists. Suppose E''/K
#    is (E')^F, the twist of E'/K by F/K. Then:
#
#    N (Cond E' x_{K} F) * D^{2} = Cond E' * Cond E''
#
#    where:
#
#    -- Cond E' is the conductor of E'/K
#    -- Cond E'' is the conductor of E''/K
#    -- D = Delta_{F/K} is the relative discriminant of F/K
#    -- E' x_{K} F is the base change of E' to F
#    -- Cond E' x_{K} F is the conductor of this base change
#    -- N = N_{F/K} is the relative norm for the extension F/K
#
#    Reference: Potential Good Reduction of Elliptic Curves,
#               Masanari Kida, J. Symbolic Computation
#               (2002) 34, 173--180.
#               doi:10.1006/jsco.2002.0555
#
#    The twist(s) E'' that we are searching for must each
#    have conductor supported over {l}, while both the norm
#    of the conductor of the base change and the conductor
#    of E are supported over {2,l}. Consequently, F/K can
#    only ramify over primes dividing 2*l. This forces F to
#    have the form F = K(\sqrt{u}), where u is a unit in
#    O_{K,S}^x, the S-units of O_K, where
#
#       S = { pp in O_K : pp | 2*l }.
#
#    Clearly u is unique only up to squares, so there are
#    only finitely many possibilities for u, and they are
#    easily enumerated: Let B = {e0, e1, ..., ek} be a
#    ZZ-basis for O_{K,S}^x. Then up to squares, the only
#    possible choices for u are elements in:
#
#       U = { \prod_{z in Z} z : Z \subseteq B, #Z != 0}.
#
#    This is a direct consequence of the m=2 case of the
#    main theorem of Kummer theory, which guarantees the
#    maximal abelian extension of K(\mu_m) is obtained by
#    adjoining mth roots of elements of K(\mu_m).
#
#    We now enumerate the elements of U, build
#    F = K(\sqrt{u}) for each u in U, and check whether
#    the curve E''/K, obtained as a quadratic twist of
#    E'/K by F/K, has good reduction away from l. If it
#    does not, we discard (E'', l) and move to the next
#    twist. If it does, we record (E''/K, l).
#
#    Notice that each step can produce only finitely many
#    objects to pass on to later steps, and so the process
#    is exhaustive and terminates in finite time. The final
#    result is a list of ``candidates'': pairs (E/K, l),
#    where
#
#    -- K is a quadratic field
#    -- E/K is a CM curve with good reduction outside l
#    -- j(E) \not\in QQ, i.e., K = Q(j(E))
#
#    This is the list HeavenlyCandidates produced below.
#
# 6. Since the condition of being heavenly at l is an
#    isogeny invariant, it is useful to sort the curves in
#    HeavenlyCandidates into isogeny classes. The final
#    part of the code distributes the candidate curves
#    into a dictionary CandidateIsogenyClass, whose
#    keys all have the form (D, l, n), where
#
#    -- D = absolute discriminant of K
#    -- l = prime at which curves could be heavenly
#    -- n = a positive integer used to distinguish
#           classes which have the same D and l
#
#    CandidateIsogenyClass[ (D, l, n) ] is a list of
#    elliptic curves. Each curve is defined over K, the
#    unique real quadratic field of absolute discriminant D,
#    and has good reduction away from l. Two curves E and E'
#    are isogenous over their common base field if and only
#    if they both appear in the same list.
#
#    One curve from each isogeny class will be passed to
#    a MAGMA script to test whether or not the curve, and
#    hence all curves from the same isogeny class, are
#    heavenly at l.
#
#
# Experimentally, default stack size is not sufficient
pari.allocatemem(10**8, 8*10**9)

# Let CM2 be the set of all d for which the quadratic number field QQ(\sqrt{d}) admits
# an elliptic curve with CM by an order with class number 2. We have

CM2 = [2, 3, 5, 6, 7, 13, 17, 21, 29, 33, 37, 41, 61, 89]

# Reference: On the number of isomorphism classes of CM elliptic curves defined
#            over a number field, Harris Daniels and Álvaro Lozano-Robledo,
#            J. Number Theory 157 (2015) 367--396.

# We create a list of discrminants from CM2:

CM2_disc = []
for D in CM2:
    if D%4 != 1:
        CM2_disc.append( 4*D )
    else:
        CM2_disc.append( D )

# We now create the 58 distinct j-invariants corresponding to elliptic curves E with CM
# by an order of class number 2, in the following form: [ mK, (j0, j1)], where
#    - K = Q(a), a is a root of the minimal polynomial mK
#    - the j-invariant is j0 + j1 * a
# j-invariants lying in the rational field QQ are not included.

jlist = []
T = polygen(ZZ, 'T')

for Delta in CM2_disc:
    mK = preferred_char_poly(Delta, T)
    K = NumberField( mK, 'a' )
    # to test rationality, we find the nontrivial element of Gal(K/Q)
    sigma = K.galois_group().gens()[0]
    for cm_datum in cm_j_invariants_and_orders(K):
        j = cm_datum[2]
        # We are only interested in j-invariants outside of QQ
        if j != sigma(j):
            j_data = [mK, j.vector()]
            jlist.append( j_data )

# For each j-invariant in jlist, we construct E/K, with K = QQ(j) and j(E) = j.

SeedCurves = []
for mK, jvec in jlist:
    K.<a> = NumberField( mK )
    j_invar = jvec[0] + jvec[1] * a
    E = EllipticCurve(j = j_invar)
    # We replace E with a global minimal model.
    E = E.global_minimal_model()
    SeedCurves.append( E )

# This concludes Step 1.

elapsed = time.time() - true_start
print(f"List of seed curves (one per j-invariant) built in {elapsed:.3f}s.", flush=True)

next_lap = time.time()

CandidatePairs = []

for E in SeedCurves:
    # We determine the CM field L, as well as the relative discriminant
    # of KL/K, then calculate the set RamQ(KL/K).
    K = E.base_ring()
    Delta_K = K.discriminant()
    mK = K.defining_polynomial()
    cm_disc = E.cm_discriminant()
    L.<b> = QuadraticField(cm_disc)
    g = L.defining_polynomial()
    KL.<omega> = K.extension(g)
    RamQ = ZZ(KL.relative_discriminant().norm()).prime_divisors()
    if len(RamQ) > 1:
        # The ramification set of KL/K is too large and
        # there is no hope of finding a twist with good
        # reduction away from a single prime. This case
        # may be skipped.
        None
    elif len(RamQ) == 1:
        # The ramification of KL/K takes place entirely
        # over the single rational prime l.
        #
        # If a twist of E is heavenly, then that twist
        # must be heavenly at l; we record (E,l)
        # for further analysis.
        l = RamQ[0]
        CandidatePairs.append( (E, l) )
    elif len(RamQ) == 0:
        # The extension KL/K is unramified (at finite places)
        # and so a priori a twist of E could be heavenly at any prime.
        # The heavenly-at-l condition, however, forces any twist
        # to have a K-rational isogeny of degree l.
        #
        # First, we find all prime numbers that are the degree
        # of a K-rational isogeny of E.
        prime_isog_degrees = []
        for phi in E.isogenies_prime_degree():
            l = phi.degree()
            if l not in prime_isog_degrees:
                prime_isog_degrees.append( l )

        # the primes in prime_isog_degrees are the only primes
        # where E could possibly be heavenly at l. We record
        # (E, l) for each l in prime_isog_degrees.
        prime_isog_degrees.sort()
        for l in prime_isog_degrees:
            CandidatePairs.append( (E, l) )

elapsed = time.time() - next_lap
print(f"Candidate pairs identified. Elapsed time {elapsed:.3f}s", flush=True)

# We have concluded Steps 2 and 3; the list CandidatePairs contains
# every pair that must be passed to Step 4.

# We replace each candidate pair (E0, l), if possible, with a pair (E1, l)
# where E1 has a minimal discriminant ideal supported over 2*l.

TwistedCandidatePairs = []

for E0, l in CandidatePairs:
    twist_exists, E1 = reduce_discriminant_by_qtwist(E0, l)
    if twist_exists:
        TwistedCandidatePairs.append( (E1, l) )

# We now check each pair (E, l) in TwistedCandidatePairs for
# quadratic twists that have good reduction away from one prime.
# There are only finitely many quadratic twists of E to check,
# described by the possible subsets of a ZZ-basis for O_{K,S}^x,
# where S = {pp in OK : pp divides 2*l}.

next_lap = time.time()

HeavenlyCandidates = []

for E, l in TwistedCandidatePairs:
    K = E.base_ring()
    S_E = [pp for pp, val in E.discriminant().factor()]
    OKSx = K.S_unit_group(S=S_E)
    twist_values = []
    B = [ K(u) for u in OKSx.gens()]
    for subcollection in Subsets(B).list():
        twist_val = prod(u for u in subcollection)
        twist_values.append( twist_val )
    for twist_val in twist_values:
        Enew = E.quadratic_twist( twist_val ).global_minimal_model()
        Snew = [p for p, val in Enew.conductor().norm().factor()]
        if len(Snew) <= 1:
            HeavenlyCandidates.append( (Enew, l) )
elapsed = time.time() - next_lap
print(f"Heavenly candidates identified. Elapsed time {elapsed:.3f}s.", flush=True)

CandidateIsogenyClass = {}
# keys for this dictionary will have the form (D, l, n), where
#     -- D is the absolute discriminant of K/Q (and this uniquely identifies K)
#     -- l is a rational prime where E could be heavenly
#     -- n is a natural number distinguishing classes where D and l are the same
# two curves in CandidateIsogenyClasses[ (D, l, n) ] will definitely be isogenous over K
# two curves in distinct lists indexed by (D, l, m) and (D, l, n) with m != n
# are definitely *not* isogenous over K.

curve_count = 0
isog_class_count = 0
for (E,l) in HeavenlyCandidates:
    curve_count += 1
    # There appear to be some pari and/or SAGE issues that cause memory errors.
    # When the used portion of the stack gets too large, we reset:
    if pari.getstack() > pari.stacksize() // 5:
        print("Clearing the PARI mechanism...", flush=True)
        pari.allocatemem(pari.stacksize(), pari.stacksizemax(), silent=True)

    D = E.base_ring().discriminant()
    n = 1
    found_isogeny_class = False
    while (D, l, n) in CandidateIsogenyClass:
        # (D,l,n) will only exist as a key if CandidateIsogenyClass[ (D,l,n) ] is nonempty
        # We grab the first curve in this isogeny class for comparison; we also time
        # the calculation.
        E0 = CandidateIsogenyClass[ (D, l, n) ][0]
        t0 = time.time()
        result = is_pair_isogenous(E, l, E0, l)
        elapsed = time.time() - t0
        if result:
            CandidateIsogenyClass[ (D, l, n) ].append(E)
            print(f" Curve #{curve_count} lies in class ({D},{l},{n}) -- isogeny took {elapsed:.3f}s to confirm.", flush=True)
            found_isogeny_class = True
            break
        n += 1
    if not found_isogeny_class:
        CandidateIsogenyClass[ (D, l, n) ] = [E]
        isog_class_count += 1
        # Also, since we just started a new isogeny class, let's report on our progress.
        print( f"Isogeny class #{str(isog_class_count).zfill(2)}, labeled ({D},{l},{n}), seeded with curve #{str(curve_count).zfill(3)}.", flush=True)

# once calculation is complete, we write a list of all isogeny classes to file
#
# heavenly-candidates-complete.sage
#
# for later use. Each line has the form
#
#    [D, l, n, mK, ec_ainv_list], where
#
#    -- (D, l, n) is the isogeny class label
#    -- K = QQ(a) = unique real quadratic field of absolute discriminant D
#    -- mK is a minimal polynomial for a (recorded in indeterminate 'T')
#    -- ec_ainv_list is a list of lists
#          -- each entry is ainv_list = [a1vec, a2vec, a3vec, a4vec, a6vec]
#             where ajvec = [aj0, aj1] is a pair of integers such that aj = aj0 + a * aj1
#             gives the K-rational a-invariants for an elliptic curve in the isogeny class.

with open('heavenly-candidates-complete.sage', 'w') as f_sage:
    for (D, l, n), curves in CandidateIsogenyClass.items():
        K = curves[0].base_ring()
        mK = K.defining_polynomial()
        ec_ainv_list = []
        for E in curves:
            ainv_list = [list(K(ai).vector()) for ai in E.a_invariants()]
            ec_ainv_list.append(ainv_list)
        f_sage.write(f"{[D, l, n, mK, ec_ainv_list]}\n")

# We also write the output to a MAGMA file, since the calculation
# to check whether E is heavenly at l will be done in MAGMA, not sage.
# That file is
#
# heavenly-candidates-complete.m
#
# in the present directory. Loading and executing this MAGMA code
# will construct an associative array (i.e., dictionary), "data",
# whose keys are lists [D, l, n] and whose value is a tuple
#
# < mK_coeffs, list_of_ainvs >
#
# where mK_coeffs will let us reconstruct the minpoly for a, and
# list_of_ainvs contains a list of lists [ [a11, a12], [a21, a22], [a31, a32], [a41, a42], [a61, a62]]
# where the a-invariants aj = aj1 + aj2 * a

with open('heavenly-candidates-complete.m', 'w') as f_magma:
    f_magma.write("data := AssociativeArray();\n\n")
    for (D, l, n), curves in CandidateIsogenyClass.items():
        K = curves[0].base_ring()
        mK = K.defining_polynomial()
        key = f"[{D}, {l}, {n}]"
        f_magma.write(f"data[{key}] := <{list(mK)}, [\n")
        for i, E in enumerate(curves):
            ainv_list = [list(K(ai).vector()) for ai in E.a_invariants()]
            cstr = ", ".join(f"[{c[0]}, {c[1]}]" for c in ainv_list)
            comma = "," if i < len(curves) - 1 else ""
            f_magma.write(f"  [{cstr}]{comma}\n")
        f_magma.write(f"]>;\n\n")
