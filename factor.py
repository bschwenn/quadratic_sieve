from sieve import *
from matrix_reducer import *
from copy import copy,deepcopy
import numpy
import math
import time
import random

TRIAL_DIVISION_CUTOFF = 100000
SMALL_PRIMES = sieve_era(1000)

def factor(n):
    return prettify(n, factor_raw(n))

def factor_raw(n, gcd_cache = set()):
    ##print("Attempting to factor {}...".format(n))
    if n < TRIAL_DIVISION_CUTOFF:
        return factor_by_division(n)

    B = find_bound(n)
    factor_base = sieve_era(B)
    smooth_squares_factor_map = {}
    pi_B = len(factor_base)
    factors = set()
    used_primes = set()

    for (x,smooth_square) in sieve_quad_poly_log(n, factor_base, B):
        # ##print("Found smooth square {} with x={}".format(smooth_square, x))
        if (smooth_square == None):
            # ##print("HEY THERE MAMA")
            factors.add(x)
            continue

        factor_map = factor_in_base_map(smooth_square, factor_base)

        smooth_squares_factor_map[x] = factor_map
        used_primes |= factor_map.keys()
        # ##print("length of smooth_squares_factor_map:", len(smooth_squares_factor_map))
        # ##print("length of  used_primes:", len(used_primes))
        if len(smooth_squares_factor_map) > len(used_primes):
            x_vector = list(smooth_squares_factor_map.keys())
            factor_matrix = []
            sorted_base = sorted(list(used_primes))

            for x in x_vector:
                factor_matrix.append([smooth_squares_factor_map[x][f] for f in sorted_base])

            start = time.time()
            left_nullspace_mat = find_dependencies(numpy.array(factor_matrix))
            # ##print("Found {} dependencies in {} seconds".format(len(left_nullspace_mat), time.time()-start))

            new_factors = check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, sorted_base, gcd_cache)

            # Find more b-smooth numbers if we didn't get any useful dependencies
            if not new_factors:
                continue

            factors |= new_factors
            #print(factors)
            result = is_fully_factored(n, factors)

            # If result returns an int it is the (probably, by Miller-Rabine) prime
            # result of dividing n by all the already found factors, i.e., its
            # the last factor
            if type(result) == int:
                factors.add(result)
                return factors
            elif result:
                return factors

            # might need to remove stuff from the map instead of just adding more
            # Update: this seems to work fine... it can just only add one dependency at a time though... we should probably add a check to avoid retesting dependencies
    return factors

def prettify(n, factors):
    result = []

    for factor in sorted(factors):
        k = 1

        power = factor
        while n % power == 0:
            if k > 1 and power in factors:
                factors.remove(power)

            k += 1
            power = factor**k

        if k == 2: # k - 1 == 1
            result.append(factor)

        else:
            power_of_prime = factor**(k-1)
            result.append("{}^{}".format(factor, k-1))

    return result

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
        #print(n, prime)
        k = 1

        while n % (prime**k) == 0:
            k += 1

        if k > 1:
            factors.add(prime**(k-1))

    #print(factors)
    return factors

def check_for_factors(n, x_vector, factor_matrix, left_nullspace_mat, factor_base, cache):
    factors = set()

    for idx,row in enumerate(left_nullspace_mat):
        sum_vec = numpy.array([ 0 for i in range(len(factor_base)) ] )

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

        # ##print("Testing equal squares x = {}, y = {}...".format(prod_of_roots, prod_of_factors))
        if prod_of_roots != prod_of_factors and prod_of_roots != ((prod_of_factors - n) % n):
            # ##print("... sanity check: x^2={}, y^2={}".format(pow(prod_of_roots,2,n),pow(prod_of_factors,2,n)))
            gcd_xy = gcd(prod_of_roots-prod_of_factors,n)
            #print("... the gcd is {}".format(gcd_xy))

            if gcd_xy == 1 or gcd_xy == n or gcd_xy in cache:
                continue

            cache.add(gcd_xy)
            #if not is_prime(gcd_xy):
            #    if gcd_xy not in cache.keys():
            #        cache[gcd_xy] = factor(gcd_xy)
            #        factors |= cache[gcd_xy]

            if is_prime(gcd_xy):
                factors.add(gcd_xy)
                #print(factors)
            else:
                #print("Recursing on {}".format(gcd_xy))
                factors |= factor_raw(gcd_xy, cache)

    return factors

def is_fully_factored(n, factors):
    prod_factors = numpy.prod(factor)

    if n == factors:
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

    #if n == 2:
    #    return True

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
        ##print("You need to input a number to factor!")
        return 1

    n = int(sys.argv[1])
    print("The factors of", n, "are", factor(n))
    return 0

if __name__ == "__main__":
   main()
