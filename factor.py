from sieve import *
from matrix_reducer import *
from copy import copy,deepcopy
import numpy
import math
import time
import random
import logging

TRIAL_DIVISION_CUTOFF = 100000
SMALL_PRIMES = sieve_era(1000)

def factor(n):
    return prettify(n, factor_raw(n))

def factor_raw(n):
    logging.info("Attempting to factor {}...".format(n))
    if n < TRIAL_DIVISION_CUTOFF:
        return factor_by_division(n)

    # If n is perfect square of primes need to know in advance, sieve won't work otherwise
    if math.floor(math.sqrt(n)) == math.sqrt(n):
        return set([int(math.sqrt(n))])

    B = find_bound(n)
    factor_base = sieve_era(B)
    smooth_squares_factor_map = {}
    pi_B = len(factor_base)
    factors = set()
    used_primes = set()

    for (x,smooth_square) in sieve_quad_poly_log(n, factor_base, B):
        logging.debug("Found smooth square {} with x={}".format(smooth_square, x))
        if (smooth_square == None):
            factors.add(x)
            continue

        factor_map = factor_in_base_map(smooth_square, factor_base)

        smooth_squares_factor_map[x] = factor_map
        used_primes |= factor_map.keys()
        if len(smooth_squares_factor_map) > len(used_primes):
            x_vector = list(smooth_squares_factor_map.keys())
            factor_matrix = []
            sorted_base = sorted(list(used_primes))

            for x in x_vector:
                factor_matrix.append([smooth_squares_factor_map[x][f] for f in sorted_base])

            start = time.time()
            left_nullspace_mat = find_dependencies(numpy.array(factor_matrix))
            logging.debug("Found {} dependencies in {} seconds".format(len(left_nullspace_mat), time.time()-start))

            for factor in check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, sorted_base):
                process_factor(factor, factors)

            result = is_fully_factored(n, factors)

            # If result returns an int it is the (probably, by Miller-Rabine) prime
            # result of dividing n by all the already found factors, i.e., its
            # the last factor
            if type(result) == int:
                factors.add(result)
                return factors
            elif result:
                return factors

    return factors

def process_factor(factor, factors):
    if factor not in factors:
        if is_prime(factor):
            factors.add(factor)
        else:
            print("Gonna factor {}".format(factor))
            factors |= factor_raw(factor)

# Format the factors in a pretty format, e.g., remove 29^2 if 29^3 is a factor
# and actually write 29^3 as 29^3 (as opposed to 24389)
def prettify(n, factors):
    result = []
    to_ignore = set()

    for factor in sorted(factors):
        if factor in to_ignore:
            continue

        logging.debug("{} is meant to be a factor of {}".format(factor, n))
        k = 1

        power = factor
        while n % power == 0:
            to_ignore.add(power)
            k += 1
            power = factor**k

        if k == 2: # k - 1 == 1
            result.append(factor)

        else:
            power_of_prime = factor**(k-1)
            result.append("{}^{}".format(factor, k-1))

    return result

# For small numbers just factor by trial division
def factor_by_division(n):
    if n == 1:
        return [1]
    if n == 2:
        return [2]
    if n == 3:
        return [3]

    primes = sieve_era(math.ceil(math.sqrt(n)))
    factors = set()

    for prime in primes:
        k = 1
        power = prime

        while n % power == 0:
            # Add all powers as factors and remove in prettify
            factors.add(power)
            k += 1
            power = prime**k

    return factors

def check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, factor_base):
    for idx,row in enumerate(left_nullspace_mat):
        sum_vec = numpy.array([ 0 for i in range(len(factor_base)) ])

        for i, val in enumerate(row):
            if val == 1:
                sum_vec += numpy.array(factor_matrix[i])

        prod_of_roots = 1
        prod_of_factors = 1
        combined_exponent_vector = [ 0 for i in range(len(factor_base)) ]

        for i in range(len(row)):
            if row[i] == 1:
                square_root = x_vector[i]
                prod_of_roots = (prod_of_roots * square_root) % n
                exponent_vector = factor_matrix[i]

                for j in range(len(exponent_vector)):
                    combined_exponent_vector[j] += exponent_vector[j]

        for j in range(len(combined_exponent_vector)):
            assert (combined_exponent_vector[j] % 2 == 0)
            prod_of_factors = (prod_of_factors * pow(factor_base[j], int(combined_exponent_vector[j]/2))) % n

        logging.debug("Testing equal squares x = {}, y = {}...".format(prod_of_roots, prod_of_factors))
        if prod_of_roots != prod_of_factors and prod_of_roots != ((prod_of_factors - n) % n):
            logging.debug("... sanity check: x^2={}, y^2={}".format(pow(prod_of_roots,2,n),pow(prod_of_factors,2,n)))
            gcd_xy = gcd(prod_of_roots-prod_of_factors,n)
            logging.debug("... the gcd is {}".format(gcd_xy))

            if gcd_xy == 1 or gcd_xy == n:
                continue

            yield gcd_xy

def is_fully_factored(n, factors):
    prod_factors = numpy.prod(list(factors))

    if n <= prod_factors:
        return True

    div_by_all_f = n

    for f in factors:
        div_by_all_f /= f

    div_by_all_f = int(div_by_all_f)

    if div_by_all_f <= 1:
        return True
    elif is_prime(div_by_all_f):
        return div_by_all_f
    else:
        return False

def is_prime(n):
    n = int(n)

    if n in SMALL_PRIMES:
        return True

    if n % 2 == 0:
        return False

    m = n - 1
    k = 0

    while m % 2 == 0:
        k += 1
        m = m / 2

    a = random.randrange(2, n-1)
    m = int(m)
    a = pow(a, m, n)

    if a == 1 or a == n - 1:
        return True

    while k > 1:
        a = pow(a,2,n)
        k -= 1

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

    logging.basicConfig(filename='factor.log', filemode='w', level=logging.DEBUG, format='%(message)s')

    n = int(sys.argv[1])
    print("The factors of", n, "are", factor(n))
    return 0

if __name__ == "__main__":
   main()
