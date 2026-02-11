/* curve_data.m */
/* Format: A list of tuples < Delta_K, ell, Coeffs > to be sent to CheckHeavenlyCurves
   Coeffs must be a sequence of 5 pairs: [ [c0, c1], ... ] 
   where [c0, c1] represents c0 + c1*a */

curves_db := [
    /* 1. The original H(47, 5).1 curve — should put some isogeny class label here */
    <5, 47, [
        [0, 0],             /* a1 */
        [0, 0],             /* a2 */
        [1, 0],             /* a3 */
        [-17578, 4136],     /* a4 = -17578 + 4136*a */
        [-962572, 324723]   /* a6 = -962572 + 324723*a */
    ]>,

    <8, 31, [
        [0, 0], 
        [0, 0], 
        [0, 0], 
        [1, 1],    /* a4 = 1 + 1*a */
        [2, -3]    /* a6 = 2 - 3*a */
    ]>,

    <13, 17, [
        [0, 0], 
        [0, 0], 
        [1, 0],    /* a3 = 1 */
        [-5, 0],   /* a4 = -5 */
        [7, 0]     /* a6 = 7 */
    ]>
];