from math import *
import sys

def selectFactorBaseLimit(n):
    # TODO - choose B
    return 1000

def primesUpTo(n):
    #TODO - implement sieve of eratosthenes
    return []

def findMSmoothSquares(n, m, base):
    x = floor(n) + 1
    smoothSquares = {}
    while x < n and len(smoothSquares) < m:
         # TODO - See if x^2 is smooth, if so store list of exponents and add to smoothSquares
         continue
    return smoothSquares

def factorInBase(m, base):
    return []

def isFactor(n, p):
    return ceil(float(n)/p) == floor(float(n)/p)

def findLinearDependenciesModTwo(smoothSquares):
    return {} # TODO -implement (maps from square to vector exponent map)

def rootsAreNotCongruentModN(n, square, factorVector):
    return False

# I think we can use a library for gcd

def factor(n, x, y):
    # return gcd(x-y, n) and n/gcd(x-y,n)
    return [0,0]

def main():
    if len(sys.argv) < 2:
        print("You need to input a number to factor!")

    n = sys.argv[1]
    # TODO - actually do things :)
    return 0

if __name__ == "__main__":
    main()
