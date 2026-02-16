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

def preferred_char_poly(Delta, indeterminate=T):
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
#       
# Experimentally, default stack size is not sufficient

pari.allocatemem()  
pari.allocatemem()

# Here is a very brief overview of the strategy. For each possible j-invariant,
# we grab /some/ integral elliptic curve E0 with that j-invariant. We insist that
# E0 is defined over the field K0 = QQ(j)
#
# We use quadratic twists to move around inside the Qbar-isomorphism class of E0,
# looking for a any representatives E that have the following properties:
#
#     - E is defined over K0
#     - E has good reduction away from one rational prime l
#
# Once we have done this for every possible j-invariant, we will have a list X of
# elliptic curves with the following property:
#
# If K1 is a quadratic field and A/K1 is a CM ellliptic curve which is heavenly at some prime l,
# and which does not have rational j-invariant, then there exists a curve E inside the list X such that
# A1 and E are K1-isogenous.
#
# We do *NOT* claim that every curve in the list X is heavenly, though we may verify this through a
# separate calculation.

# Let CM2 be the set of all d for which the quadratic number field QQ(\sqrt{d}) admits
# an elliptic curve with CM by an order with class number 2. We have

CM2 = [2, 3, 5, 6, 7, 13, 17, 21, 29, 33, 37, 41, 61, 89]

# Reference: 
#
# On the number of isomorphism classes of CM elliptic curves defined over a number field, 
# by Harris Daniels and Álvaro Lozano-Robledo, J. Number Theory 157 (2015) 367--396.

# we create a list of discrminants from CM2:

CM2_disc = []
for D in CM2:
    if D%4 != 1:
        CM2_disc.append( 4*D )
    else:
        CM2_disc.append( D )

# We now create the 58 distinct j-invariants corresponding to elliptic curves E with CM
# by an order of class number 2, in the following form: [ mK, (j0, j1)], where
#    - K has minimal polynomial mK
#    - if a is a root of K, then the j-invariant is j0 + j1 * a
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
        if j != sigma(j):
            # now sure j is not rational
            j_data = [mK, j.vector()]
            jlist.append( j_data )

def cm_disc_to_field_disc_and_conductor(cm_disc, class_no = 2):
    r'''
    Every CM order has the form O = ZZ + f * O_L, and is determined by its CM discriminant, cm_disc.
    
    INPUT:   
        - cm_disc, the cm discriminant of some CM order O.
        - class_no, the specified class number for the order
    OUTPUT:
        - [D, f], where
        - L = QuadraticField(D) is the CM field containing O,
        - f is the CM conductor such that O = ZZ + f * O_L.'''
    
    # The built-in function cm_orders() contains the relevant information; we just search
    # for the correct values.

    for D, f in cm_orders(class_no):
        if D * f ** 2 == cm_disc:
            return [D, f]