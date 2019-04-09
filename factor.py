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
MILLER_RABIN_NUM_ROUNDS = 3

def factor(n):
    return prettify(n, factor_raw(n))

# Obtain nontrivial factors. Result may contain, for example, 29^2 and 29^3 if both
# divide n. This is accounted for in prettification.
def factor_raw(n):
    logging.info("Attempting to factor {}...".format(n))

    # Don't try to factor prime numbers
    if is_prime(n):
        return [n]

    # If n is small, we can factor quickly with trial division
    if n < TRIAL_DIVISION_CUTOFF:
        factors = factor_by_division(n)
        logging.debug("Found factors {} for {} by trial division".format(factors, n))
        return factor_by_division(n)

    # Get factors less than log(n)
    small_factors = get_small_prime_factors(n)

    if small_factors:
        result = is_fully_factored(n, small_factors)

        if type(result) == int:
            small_factors.add(result)

        if result:
            return small_factors
        else:
            reduced_n = n

            for f in small_factors:
                reduced_n = reduced_n//f

            # Recurse after dividing out small factors
            large_factors = factor_raw(reduced_n)
            return small_factors | large_factors

    # If n is perfect square of primes need to know in advance
    if math.floor(math.sqrt(n)) == math.sqrt(n) and is_prime(int(math.sqrt(n))):
        return set([int(math.sqrt(n))])

    # Setup sieving
    B = find_bound(n)
    factor_base = sieve_era(B)
    # This map tracks the value of x in the polynomial sequence that led to the
    # smooth square, and the associated factorization of the smooth square.
    smooth_squares_factor_map = {}
    # Keep track of the primes in the factor base that we actually use in some
    # factorization, so as to reduce the size of the matrix during the Gaussian
    # elimination phase of the algorithm.
    used_primes = set()
    # Prime (or power of prime, for now) factors
    factors = set()

    # Obtain smooth primes one at a time from the sieving function
    for (x,smooth_square) in sieve_quad_poly_log(n, factor_base, B):
        # If we luck upon a factor in sieving smooth_square will be None
        if (smooth_square == None):
            factors.add(x)
            continue

        logging.debug("Found smooth square {} with x={}".format(smooth_square, x))
        # Find the factorization of the smooth square in the factor base.
        # Experimentally, this was faster than having the sieving function
        # Maintain the factorizations as it determined smooth squares.
        factor_map = factor_in_base_map(smooth_square, factor_base)

        # factor_map maps the prime in the factor base to its exponent
        # in the factorization.
        smooth_squares_factor_map[x] = factor_map
        used_primes |= factor_map.keys()

        logging.debug("Used {} primes in {} smooth squares".format(len(used_primes), len(smooth_squares_factor_map)))
        # Once we more smooth squares than primes used in their factorizations,
        # we are guaranteed to find a dependency.
        if len(smooth_squares_factor_map) > len(used_primes):
            x_vector = list(smooth_squares_factor_map.keys())
            factor_matrix = []
            sorted_base = sorted(list(used_primes))

            # Construct the factor matrix (entries being the exponents in the factorization for the
            # corresponding prime in the factor base for each x).
            for x in x_vector:
                factor_matrix.append([smooth_squares_factor_map[x][f] for f in sorted_base])

            start = time.time()
            # Find the dependenc(y)(ies) (guaranteed to exist)
            left_nullspace_mat = find_dependencies(numpy.array(factor_matrix))
            logging.debug("Found {} dependencies in {} seconds".format(len(left_nullspace_mat), time.time()-start))

            # Obtain nontrivial factors from dependecies where possible
            for factor in check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, sorted_base):
                process_factor(n, factor, factors)

            # Check if we have fully factored n
            result = is_fully_factored(n, factors)

            # If is_fully_factored returns an int it is the (probably, by Miller-Rabine) prime result
            # of dividing n by all the already found factors, i.e., its the last factor.
            if type(result) == int:
                factors.add(result)
                return factors
            elif result:
                return factors

    return factors

# Check if the factor is new, and if it can be factored further, and
# then handle appropriately
def process_factor(n, factor, factors):
    if factor not in factors:
        if is_prime(factor):
            logging.debug("{} is apparently prime".format(factor))
            factors.add(find_highest_power_factor(n, factor))
        else:
            factors |= factor_raw(factor)

# Format the factors in a pretty format, e.g., remove 29^2 if 29^3 is a factor
# and actually write 29^3 as 29^3 (as opposed to 24389).
def prettify(n, factors):
    # List of factors in nice format
    result = []
    # Ignore powers of primes that occur later in factors, since we check for
    # the highest power of primes earlier in the list which divide n. (This
    # is why we sort the factors.)
    to_ignore = set()

    for factor in sorted(factors):
        if factor in to_ignore:
            continue

        logging.debug("{} is meant to be a factor of {}".format(factor, n))
        k = 1
        power = factor

        # Find highest dividing power
        while n % power == 0:
            to_ignore.add(power)
            k += 1
            power = factor**k

        if k == 2: # highest dividing power is to the first
            result.append(factor)
        else:
            result.append("{}^{}".format(factor, k-1))

    return result

# Find the small prime factors of n (<= log(n)) by trivial division.
def get_small_prime_factors(n):
    logn = int(math.log(n))
    return get_factors_up_to(n, logn)

# Find all prime factors of n up to limit.
def get_factors_up_to(n, limit):
    primes = sieve_era(limit)
    factors = set()

    for prime in primes:
        power = find_highest_power_factor(n, prime)

        if power:
            factors.add(power)

    return factors

# Find the highest power of a prime dividing n (just the power, not exponent).
def find_highest_power_factor(n, prime):
    k = 1
    power = prime
    old_power = None

    while n % power == 0:
        k += 1
        old_power = power
        power = prime**k

    return old_power

# Factor by trial division, used for small n.
def factor_by_division(n):
    if n == 1:
        return [1]
    if n == 2:
        return [2]

    return get_factors_up_to(n, n)

# Given matrix of dependencies (left_nullspace_mat), check if these dependencies
# yield nontrivial factors of n.
def check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, factor_base):
    for idx,row in enumerate(left_nullspace_mat):
        # The product the square roots (the x's in the polynomial sequence)
        prod_of_roots = 1
        # The product of the elements in the factorizations (raised to
        # the appropriate powers)
        prod_of_factors = 1

        # Calculate the appropriate exponents
        combined_exponent_vector = [ 0 for i in range(len(factor_base)) ]

        for i in range(len(row)):
            if row[i] == 1:
                square_root = x_vector[i]
                prod_of_roots = (prod_of_roots * square_root) % n
                exponent_vector = factor_matrix[i]

                for j in range(len(exponent_vector)):
                    combined_exponent_vector[j] += exponent_vector[j]

        # Calculate the product of the factors raised to appropriate powers
        for j in range(len(combined_exponent_vector)):
            assert (combined_exponent_vector[j] % 2 == 0)
            prod_of_factors = (prod_of_factors * pow(factor_base[j], int(combined_exponent_vector[j]/2))) % n

        # Check if the x is not congruent to y in x^2 cong y^2 mod n
        logging.debug("Testing equal squares x = {}, y = {}...".format(prod_of_roots, prod_of_factors))
        if prod_of_roots != prod_of_factors and prod_of_roots != ((prod_of_factors - n) % n):
            logging.debug("... sanity check: x^2={}, y^2={}".format(pow(prod_of_roots,2,n),pow(prod_of_factors,2,n)))
            gcd_xy = gcd(prod_of_roots-prod_of_factors,n)
            logging.debug("... the gcd is {}".format(gcd_xy))

            if gcd_xy == 1 or gcd_xy == n:
                continue

            # ... if so, we found a nontrivial factor
            yield gcd_xy

# Check if the factors supplied lead to a complete factorization of n
def is_fully_factored(n, factors):
    prod_factors = 1

    for factor in factors:
        prod_factors = (prod_factors * factor) % n

    if prod_factors == 0:
        return True

    leftover = n / prod_factors

    if leftover == int(leftover) and is_prime(int(leftover)):
        return int(leftover)
    else:
        return False

# Miller-Rabin primality test
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

# Command-line interface stuff
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
