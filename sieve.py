import math
import sys
from collections import defaultdict

def sieve_era(B): # returns primes <= b
    l = [True] * (B+1)
    l[0] = l[1] = False
    for i in range(2,int(B**0.5)+1):
        if l[i]:
            for a in range(i*i,B+1,i):
                l[a] = False
    return [i for i,j in enumerate(l) if j]

def pi(B):
    return len(sieve_era(B))-1

def get_sol(n,p): #gets solutions to x^2-n = 0 mod p
    sol = []
    for i in range(p):
        if (((i*i)-n) % p == 0):
            sol.append(i)
    return sol

def sieve_quad_poly_log(n, p_set, B): #n should be odd
    x_start = int(math.ceil(n**0.5))+10  # bounds of sieve
    top_bound = int(math.ceil(n**0.5))+16*B

    orig = [i*i-n for i in range(x_start, top_bound)]  # preserved list
    l = [math.log(i*i - n) for i in range(x_start,top_bound)]  # list of all the corresponding x^2-n that will be modified

    for p in p_set:
        sols = get_sol(n,p) # solutions x for x^2-n=0 mod p

        if p == 2: # solutions to x^2-n = 0 mod 2: if n is even, x = 0 mod 2, if n is odd, x = 1 mod 2
            index = 0
            for x in range(1+2*math.ceil((x_start-1)/2), top_bound, 2):
                k = 1

                # Account for the number being divisible by the prime more than once
                while orig[x-x_start] % 2**k == 0:
                    l[x-x_start] -= math.log(2)
                    k += 1

                if (l[x-x_start]<0.1):
                    yield (x, orig[x-x_start])

        elif len(sols) == 1:  # p | n
            yield (p, None)
        elif len(sols) == 2:  # two solutions, sieve
            for x in range(sols[0]+p*math.ceil((x_start-sols[0])/p), top_bound, p):
                num = x-x_start
                k = 1

                # Account for the number being divisible by the prime more than once
                while orig[x-x_start] % p**k == 0:
                    l[num] -= math.log(p)
                    k += 1

                if l[num] < 0.1:
                    yield (x, orig[x-x_start])

            for x in range(sols[1]+p*math.ceil((x_start-sols[1])/p),top_bound,p):
                num = x-x_start
                k = 1

                # Account for the number being divisible by the prime more than once
                while orig[num] % p**k == 0:
                    l[num] -= math.log(p)
                    k += 1

                if l[num] < 0.1:
                    yield (x, orig[x-x_start])

def factor_in_base(n,p_set): # return exponent vector for n, w/ primes up to B through trial division
    ans = []

    for i in p_set:
        count = 0

        while n % i == 0:
            n /= i
            count += 1

        ans.append(count)

    return ans

def factor_in_base_map(n, p_set):
    ans = defaultdict(lambda: 0)

    for i in p_set:
        count = 0
        exp = 1

        while n % i**exp == 0:
            exp += 1

        if exp-1 != 0:
            ans[i] = exp-1
    return ans

def gcd(a, b):
    return gcd(b, a % b) if a % b else b

def find_bound(n):
    return int(math.exp((2**(-0.5)*(math.log(2*n)*math.log(math.log(2*n)))**0.5)))
