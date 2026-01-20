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

# We now collect data to represent every possible isomorphism 

EllipticCurveCM order of class number 2, , with the following structure:
#   - char_poly is the minimal polynomial for the field of definition K = QQ(a),
#       chosen in a standard way and such that char_poly(a) == 0.
#   - 
#   - Delta_L is the absolute discriminant of the CM field L
#   - f is the CM conductor of the desired order (e.g., O = ZZ + f * O_L)
#   - j 




def preferred_char_poly(Delta, indeterminate=T):
    r'''
    For a positive real integer Delta, returns a standardized form
    for the minimal polynomial of the unique quadratic field K of absolute
    discriminant Delta, in the specified indeterminate.'''

    if Delta % 4 == 0:
        return indeterminate^2 - (Delta/4)
    elif:
        return indeterminate^2 - indeterminate - (Delta - 1)/4
    else:
        raise Exception("Delta is not the discriminant of a quadratic field.")


def cm_disc_to_field_disc_and_conductor(cm_disc, class_no = 2):
    r'''
    Every CM order has the form O = ZZ + f * O_L, and is determined by its CM discriminant, cm_disc.
    
    INPUT:   
        - cm_disc, the cm discriminant of some CM order O.
    OUTPUT:
        - [D, f], where
        - L = QuadraticField(D) is the CM field containing O,
        - f is the CM conductor such that O = ZZ + f * O_L.'''
    
    # The built-in function cm_orders() contains the relevant information; we just search
    # for the correct values.

    for D, f in cm_orders(class_no):
        if D * f ** 2 == cm_disc:
            return [D, f]