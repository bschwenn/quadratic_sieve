from sieve import *
from matrix_reducer import *
import numpy
import math

def factor(n):
    B = find_bound(n)
    factor_base = sieve_era(B)
    smooth_squares_factor_map = {}
    pi_B = len(factor_base)
    factors = set()

    for (x,smooth_square) in sieve_quad_poly_log(n, factor_base, B):
        print(x, smooth_square)
        if (smooth_square == None):
            factors.add(x)
            continue
#        print("Square:" + str(smooth_square))
        smooth_squares_factor_map[x] = factor_in_base(smooth_square, factor_base)
        #print(str(smooth_square) + " " + str(smooth_squares_factor_map[x]))

        if len(smooth_squares_factor_map) > pi_B:
            # print([v for k,v in smooth_squares_factor_map.items()])
            x = [v for k, v in smooth_squares_factor_map.items()]
            left_nullspace_mat = find_dependencies(numpy.array(x))
            print(left_nullspace_mat)
            # print(np.matmul(left_nullspace_mat, numpy.array(x)))

            new_factors = check_for_factors(n, smooth_squares_factor_map, left_nullspace_mat, factor_base)
            print("New factors: ", new_factors)
            factors |= new_factors

            result = is_fully_factored(n, factors)

            if type(result) == int:
                factors.add(result)
                return factors
            elif result:
                return factors

            # might need to remove stuff from the map instead of just adding more
            print("Factors: ", factors)

    print(len(smooth_squares_factor_map), " ", pi_B)
    return factors

def check_for_factors(n, smooth_squares_factor_map, left_nullspace_mat, factor_base):
    roots_of_squares = list(smooth_squares_factor_map.keys())
    factors = set()
    for row in left_nullspace_mat:
        prod_of_roots = 1
        prod_of_factors = 1

        for i in range(len(row)):
            if row[i] == 1:
                square_root = roots_of_squares[i]
                prod_of_roots *= square_root
                exponent_vector = smooth_squares_factor_map[square_root]

                for j in range(len(exponent_vector)):
                    prod_of_factors *= pow(factor_base[j], exponent_vector[j], n)
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


def is_prime(n):
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
    print(factor(n))
    return 0

if __name__ == "__main__":
   main()

