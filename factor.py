from sieve import *
from matrix_reducer import *
from copy import copy,deepcopy
import numpy
import math

def factor(n):
    B = find_bound(n)
    factor_base = sieve_era(B)
    smooth_squares_factor_map = {}
    pi_B = len(factor_base)
    factors = set()

    for (x,smooth_square) in sieve_quad_poly_log(n, factor_base, B):
        if (smooth_square == None):
            factors.add(x)
            continue
        smooth_squares_factor_map[x] = factor_in_base(smooth_square, factor_base)

        if len(smooth_squares_factor_map) > pi_B:
            x_vector = list(smooth_squares_factor_map.keys())
            factor_matrix = []

            for x in x_vector:
                factor_matrix.append(smooth_squares_factor_map[x])


            print(numpy.array(factor_matrix))
            left_nullspace_mat = find_dependencies(numpy.array(factor_matrix))
            print(left_nullspace_mat)
            new_factors = check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, factor_base)
            factors |= new_factors

            result = is_fully_factored(n, factors)

            if type(result) == int:
                factors.add(result)
                return factors
            elif result:
                return factors

            # might need to remove stuff from the map instead of just adding more

    print(len(smooth_squares_factor_map), " ", pi_B)
    return factors

def check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, factor_base):
    factors = set()

    print(numpy.matmul(left_nullspace_mat, factor_matrix))
    for idx,row in enumerate(left_nullspace_mat):
        prod_of_roots = 1
        prod_of_factors = 1
        combined_exponent_vector = [ 0 for i in range(len(factor_base)) ]
        print(left_nullspace_mat)
        print(x_vector)
        for i in range(len(row)):
            if row[i] == 1:
                square_root = x_vector[i]
                prod_of_roots *= square_root
                exponent_vector = factor_matrix[i]

                for j in range(len(exponent_vector)):
                    #prod_of_factors *= pow(factor_base[j], exponent_vector[j], n)
                    if j == 5:
                        print("Here: ", exponent_vector[j], i, j, idx)
                    combined_exponent_vector[j] += exponent_vector[j]

        for j in range(len(combined_exponent_vector)):
            print(combined_exponent_vector[j])
            # assert (combined_exponent_vector[j] % 2 == 0)
            prod_of_factors *= pow(factor_base[j], int(combined_exponent_vector[j]/2), n)

        prod_roots_residue = prod_of_roots % n
        prod_factors_residue = int(prod_of_factors**0.5) % n
        print("x: ", prod_roots_residue, " y: ", prod_factors_residue)
        if prod_roots_residue != prod_factors_residue and prod_roots_residue != (prod_factors_residue - n):
            factors.add(gcd(prod_roots_residue-prod_factors_residue, n))


    return factors

def is_fully_factored(n, factors):
    prod_factors = numpy.prod(factor)

    if n == factors:
        return True

    div_by_all_f = n

    for f in factors:
        div_by_all_f /= f

    div_by_all_f = int(div_by_all_f)

    if is_prime(div_by_all_f):
        print(div_by_all_f)
        return div_by_all_f
    else:
        return False

def is_prime(m):
    return True # TODO - Miller Rabine

def main():
    if len(sys.argv) < 2:
        print("You need to input a number to factor!")
        return 1

    n = int(sys.argv[1])
    print(factor(n))
    return 0
if __name__ == "__main__":
   main()

