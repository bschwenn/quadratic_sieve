from sieve import *
from matrix_reducer import *

def factor(n):
    B = find_bound(n)
    factor_base = sieve_era(B)
    smooth_squares_factor_map = {}
    pi_B = len(factor_base)
    factors = set()

    for square in sieve_quad_poly_log(n, factor_base, B):
        smooth_squares_factor_map[square] = factor_in_base(square, factor_base)

        if len(smooth_squares_factor_map) > pi_B:
            left_nullspace_mat = find_dependencies([v for k,v in smooth_squares_factor_map])
            factors |= check_for_factors(n, smooth_squares_factor_map, left_nullspace_mat, factor_base)

            if is_fully_factored(n, factors):
                return factors

            # might need to remove stuff from the map instead of just adding more


def check_for_factors(n, smooth_squares_factor_map, left_nullspace_mat, factor_base):
    roots_of_squares = smooth_squares_factor_map.keys()
    factors = set()

    for row in left_nullspace_mat:
        prod_of_roots = 1
        prod_of_factors = 1

        for i in range(len(row)):
            if row[i] == 1:
                square_root = roots_of_squares[i]
                prod_of_roots *= square_root
                exponent_vector = smooth_squares_factor_map[square_roots]

                for j in len(exponent_vector):
                    prod_of_factors *= (factor_base[j]**(exponent_vector[j]/2))

        prod_roots_residue = prod_of_roots % n
        prod_factors_residue = prod_of_factors % n

        if prod_roots_residue != prod_factors_residue and prod_roots_reside != (prod_factors_residue - n):
            factors.add(gcd(prod_roots_residue-prod_factors_residue, n))

    return factors

def is_fully_factored(n, factors):
    prod_factors = reduce((lambda x, y: x * y), factor)
    return n == factors
