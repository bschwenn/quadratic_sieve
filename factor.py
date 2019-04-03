from sieve import *
from matrix_reducer import *
from copy import copy,deepcopy
import numpy
import math
import time

def factor(n):
    print("Attempting to factor {}...".format(n))
    B = find_bound(n)
    factor_base = sieve_era(B)
    smooth_squares_factor_map = {}
    pi_B = len(factor_base)
    factors = set()
    used_primes = set()

    for (x,smooth_square) in sieve_quad_poly_log(n, factor_base, B):
        #print("Found smooth square {} with x={}".format(smooth_square, x))
        if (smooth_square == None):
            factors.add(x)
            continue

        factor_map = factor_in_base_map(smooth_square, factor_base)
        smooth_squares_factor_map[x] = factor_map # factor_in_base(smooth_square, factor_base)
        used_primes |= factor_map.keys()

        if len(smooth_squares_factor_map) > len(used_primes): # pi_B:
            x_vector = list(smooth_squares_factor_map.keys())
            factor_matrix = []
            sorted_base = sorted(list(used_primes))

            for x in x_vector:
                factor_matrix.append([smooth_squares_factor_map[x][f] for f in sorted_base])

            start = time.time()
            left_nullspace_mat = find_dependencies(numpy.array(factor_matrix))
            print("Found {} dependencies in {} seconds".format(len(left_nullspace_mat), time.time()-start))

            new_factors = check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, sorted_base) # factor_base)
            if not new_factors:
                continue

            factors |= new_factors

            result = is_fully_factored(n, factors)

            if type(result) == int:
                factors.add(result)
                return factors
            elif result:
                return factors

            # might need to remove stuff from the map instead of just adding more

    return factors

def check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, factor_base):
    factors = set()

    #print("Product")
    for idx,row in enumerate(left_nullspace_mat):
        #print("Dependency")
        #print(row)
        sum_vec = numpy.array([ 0 for i in range(len(factor_base)) ] )
        for i, val in enumerate(row):
            if val == 1:
                #print(sum_vec)
                sum_vec += numpy.array(factor_matrix[i])
        #print("Sum vec:", sum_vec)
        prod_of_roots = 1
        prod_of_factors = 1
        combined_exponent_vector = [ 0 for i in range(len(factor_base)) ]
        #print(left_nullspace_mat)
        #print(x_vector)
        for i in range(len(row)):
            if row[i] == 1:
                square_root = x_vector[i]
                #print("Using x_{}={} (as in the polynomial sequence)".format(i, square_root))
                prod_of_roots = (prod_of_roots * square_root) % n
                exponent_vector = factor_matrix[i]

                for j in range(len(exponent_vector)):
                    #prod_of_factors *= pow(factor_base[j], exponent_vector[j], n)
                    #if j == 5:
                    #    print("Here: ", exponent_vector[j], i, j, idx)
                    combined_exponent_vector[j] += exponent_vector[j]
                    #print("x_i has factor {}^{}".format(factor_base[j], exponent_vector[j]))

        for j in range(len(combined_exponent_vector)):
            assert (combined_exponent_vector[j] % 2 == 0)
            prod_of_factors = (prod_of_factors * pow(factor_base[j], int(combined_exponent_vector[j]/2))) % n

        prod_roots_residue = prod_of_roots
        prod_factors_residue = prod_of_factors
        print("Testing equal squares x = {}, y = {}...".format(prod_roots_residue, prod_factors_residue))
        if prod_roots_residue != prod_factors_residue and prod_roots_residue != (prod_factors_residue - n):
            print("... sanity check: x^2={}, y^2={}".format(pow(prod_of_roots,2,n),pow(prod_of_factors,2,n)))
            gcd_xy = gcd(prod_roots_residue-prod_factors_residue,n)
            print("... the gcd is {}".format(gcd_xy))

            if gcd_xy == 1 or gcd_xy == n:
                continue

            factors.add(gcd_xy)

    return factors

def is_fully_factored(n, factors):
    prod_factors = numpy.prod(factor)

    if n == factors:
        return True

    div_by_all_f = n

    for f in factors:
        div_by_all_f /= f

    div_by_all_f = int(div_by_all_f)

    if div_by_all_f == 1:
        return True
    elif is_prime(div_by_all_f):
        return div_by_all_f
    else:
        return False

def is_prime(n):
    return True
    m = n - 1
    k = 0
    while m % 2 == 0:
        k += 1
        m = m / 2
    a = 2
    a = (a ** m) % n
    if a == 1 or a == n - 1:
        return True
    while k > 1:
        a = (a ** 2) % n
        if a == 1:
            return False
        if a == n - 1:
            return True
    if a == n - 1:
        return True
    return False

def main():
    if len(sys.argv) < 2:
        print("You need to input a number to factor!")
        return 1

    n = int(sys.argv[1])
    print("The factors of", n, "are", factor(n))
    return 0
if __name__ == "__main__":
   main()

