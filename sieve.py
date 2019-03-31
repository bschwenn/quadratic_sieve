from math import *
import sys

def selectFactorBaseLimit(n):
    # TODO - choose B
    return 1000

def primesUpTo(n):
    #TODO - implement sieve of eratosthenes
    return []
def sieve_era(B): # returns primes <= b
    l = [True] * (B+1) 
    l[0] = l[1] = False
    for i in range(2,int(B**0.5)+1):
        if l[i]:
            for a in range(i*i,B+1,i):
                l[a] = False
    return [i for i,j in enumerate(l) if j]

def sieve_quad(n,B): #return B-smooth numbers up to n, with the exponent factorizations 
    l = list(range(n+1))
    p_set = sieve_era(B) #(primes <=5)
    expfac = [ [0]*len(p_set) for i in range(n+1)]

    for i in range(2, n+1):
        if l[i] in p_set:
            index = p_set.index(l[i])
            prime = l[i]
            for j in range(i,n+1,i):
                num = j
                while (l[j] % prime == 0):
                    expfac[num][index] += 1
                    l[j] /= prime
                    
    return [(i,expfac[i]) for i,j in enumerate(l) if j==1]

def pi(B):
    return len(sieve_era(B))-1

def get_sol(x,n,p): #gets solutions to x^2-n = 0 mod p
    sol = []
    for i in range(p):
        if (i*i-n % p == 0):
            sol.append(i)
    return sol

def exp(n,B): #return exponent vector for n, w/ primes up to B through trial division
    p_set = sieve_era(B)
    ans = []
    for i in p_set:
        count = 0
        while (n % i == 0):
            n /= i
            count += 1
        ans.append(count)
    return ans


def findMSmoothSquares(n, m, base):
    x = floor(sqrt(n)) + 1
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
def gcd(a,b):
    return gcd(b,a%b) if a%b else b

def factor(n, x, y):
    # return gcd(x-y, n) and n/gcd(x-y,n)
    return [0,0]

def main():
    if len(sys.argv) < 2:
        print("You need to input a number to factor!")
        return 1

    n = sys.argv[1]
    # TODO - actually do things :)
    return 0

if __name__ == "__main__":
    main()
