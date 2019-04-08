
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

'''
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
'''

def pi(B):
    return len(sieve_era(B))-1

# def get_sol(n,p): #gets solutions to x^2-n = 0 mod p
#     sol = []
#     for i in range(p):
#         if (((i*i)-n) % p == 0):
#             sol.append(i)
#     # print(sol)
#     return sol

def has_sq_root(a, p):
    """
    Computes a^((p-1)/2) and returns 1 or -1 depending on if a has a square root mod p
    """
    x  = pow(a, (p - 1) // 2, p)
    return -1 if x == p - 1 else 1


def get_sol(n, p):
    # Algo got from wikipedia, just translated their algo into code
    # First deal w special cases
    if has_sq_root(n, p) != 1:
        return []
    elif n == 0:
        return []
    elif p == 2:
        return []
    elif p % 4 == 3:
        soln = pow(n, (p + 1) // 4, p)
        return [soln, p-soln] 

    # Now the real algorithm
    # first write p-1 as q*2^s with q odd
    m, q = 0, p-1
    while(q % 2 == 0):
        q//=2
        m+=1
    z = 2

    # searching for z which is a quadratic non-residue
    while has_sq_root(z, p) != -1:
        z += 1    

    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q+1)//2, p)

    while True:
        i = 0
        while i < m-1:
            if t == 1:
                break
            t = pow(t, 2, p)
            i += 1
        if i == 0:
            return [r, p-r]

        b = pow(c, pow(2, m-i-1, p-1), p)
        m = i
        c = pow(b, 2, p)
        t = (t*c) % p 
        r = (r*b) % p 




'''
def sieve_quad_poly(n,B): #n should be odd
    """examine numbers x^2-n for B-Smooth values, x runs integers from ceil(sqrt(n))
    returns tuple of integer, exponent vector if B-smooth in sequence x^2-n, running from x_start to top_bound
    examine numbers x^2-n for B-Smooth values, x runs integers from ceil(sqrt(n))
    returns tuple of integer, exponent vector if B-smooth in sequence x^2-n, running from x_start to top_bound"""

    x_start = int(math.ceil(n**0.5)) #bounds of sieve
    top_bound = 2 * x_start #bounds of sieve
    #print(x_start)
    print(top_bound)
    p_set = sieve_era(B)
    expfac=[ [0]*len(p_set) for i in range(x_start, top_bound)] #exponent vector

    orig = [i*i-n for i in range(x_start, top_bound)] #preserved list
    l = [i*i-n for i in range(x_start,top_bound)] #list of all the corresponding x^2-n that will be modified

    for p in p_set:
        sols = get_sol(n,p) #solutions x for x^2-n=0 mod p

        if (p == 2): #solutions to x^2-n = 0 mod 2: if n is even, x = 0 mod 2, if n is odd, x = 1 mod 2
            index = 0
            for j in range(1+2*math.ceil((x_start-1)/2), top_bound, 2):
                num = j-x_start
                while(l[j-x_start] % 2 == 0):
                    expfac[num][index]+=1
                    l[j-x_start] /= 2
        elif len(sols)==1: #p | n
            return p
        elif len(sols)==2: #two solutions, sieve
            index = p_set.index(p)

            for j in range(sols[0]+p*math.ceil((x_start-sols[0])/p),top_bound,p):
                num = j-x_start
                while (l[j-x_start] % p == 0):
                    l[j-x_start] /= p
                    expfac[num][index] += 1
            for j in range(sols[1]+p*math.ceil((x_start-sols[1])/p),top_bound,p):
                num = j-x_start
                while(l[j-x_start] % p == 0):
                    l[j-x_start] /= p
                    expfac[num][index] += 1

    return [(orig[i],expfac[i]) for i,j in enumerate(l) if j == 1]
'''

def sieve_quad_poly_log(n, p_set, B): #n should be odd
    x_start = int(math.ceil(n**0.5))+1  # bounds of sieve
#    top_bound = math.ceil(1.005 * x_start)  # bounds of sieve
    top_bound = int(math.ceil(n**0.5))+16*B

    # p_set = sieve_era(B)
    orig = [i*i-n for i in range(x_start, top_bound)]  # preserved list
    l = [math.log(i*i - n) for i in range(x_start,top_bound)]  # list of all the corresponding x^2-n that will be modified

#    b_smooth_list = []
    for p in p_set:
        sols = get_sol(n,p) # solutions x for x^2-n=0 mod p

        if p == 2: # solutions to x^2-n = 0 mod 2: if n is even, x = 0 mod 2, if n is odd, x = 1 mod 2
            index = 0
            for x in range(1+2*math.ceil((x_start-1)/2), top_bound, 2):
                k = 1

                while orig[x-x_start] % 2**k == 0:
                    l[x-x_start] -= math.log(2)
                    k += 1

                if (l[x-x_start]<0.1):
                    # b_smooth_list.append(oreig[j-x_start])
                    #print(str(x) + " " + str(orig[x-x_start]))
                    yield (x, orig[x-x_start])
                    # print(orig[j-x_start])
                #if len(b_smooth_list)>pi(B):
                    # return b_smooth_list

        elif len(sols) == 1:  # p | n
            #print(str(p) + " " + str(orig[x-x_start]))
            #print(sols, n)
            yield (p, None)
        elif len(sols) == 2:  # two solutions, sieve

            for x in range(sols[0]+p*math.ceil((x_start-sols[0])/p), top_bound, p):
                num = x-x_start
                k = 1
                while orig[x-x_start] % p**k == 0:
                    l[num] -= math.log(p)
                    k += 1
                if l[num] < 0.1:
                    #b_smooth_list.append(orig[num])
                    #print(str(x) + " " + str(orig[x-x_start]))
                    yield (x, orig[x-x_start])
                    #print(orig[num])
                #if len(b_smooth_list) > pi(B):
                #    return b_smooth_list

            for x in range(sols[1]+p*math.ceil((x_start-sols[1])/p),top_bound,p):
                num = x-x_start
                k = 1
                while orig[num] % p**k == 0:
                    l[num] -= math.log(p)
                    k += 1
                if l[num] < 0.1:
                    #b_smooth_list.append(orig[num])
                    yield (x, orig[x-x_start])
                    #print(orig[num])
                #if len(b_smooth_list) > pi(B):
                #    return b_smooth_list

    # return b_smooth_list


def factor_in_base(n,p_set): # return exponent vector for n, w/ primes up to B through trial division
    # p_set = sieve_era(B)
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
        while n % i**exp == 0: # n % i == 0:
            #n /= i
            #count += 1
            exp += 1
        if exp-1 != 0:
            ans[i] = exp-1
    return ans

def roots_congruent_mod_n(n, square, p_set, exponent_vector):
    pass


def gcd(a, b):
    return gcd(b, a % b) if a % b else b

def factor(n, x, y):
    # TODO - return gcd(x-y, n) and n/gcd(x-y,n)
    return [0, 0]


def find_bound(n):
    return int(math.exp((2**(-0.5)*(math.log(2*n)*math.log(math.log(2*n)))**0.5)))
