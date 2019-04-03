from sieve import *
import unittest

class SieveTestCase(unittest.TestCase):

    def test_factor_in_base_map_small(self):
        # large m
        m = 97157097446703
        p_set = sieve_era(1580) # will give all the primes in its factorization
        factor_map = factor_in_base_map(m, p_set)

        prod = 1
        for k,v in factor_map.items():
            prod *= k**v

        self.assertEqual(m, prod)

    def test_factor_in_base_map_large(self):
        # large m
        m = 916601041922090486400
        p_set = sieve_era(2684) # will give all the primes in its factorization
        factor_map = factor_in_base_map(m, p_set)

        prod = 1
        for k,v in factor_map.items():
            prod *= k**v

        self.assertEqual(m, prod)
