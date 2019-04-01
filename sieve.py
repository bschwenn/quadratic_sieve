import math
import sys


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


def get_sol(n,p): #gets solutions to x^2-n = 0 mod p
    sol = []
    for i in range(p):
        if (i*i-n % p == 0):
            sol.append(i)
    return sol


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

def sieve_quad_poly_log(n,p_set, B): #n should be odd

    x_start = int(math.ceil(n**0.5))  # bounds of sieve
#    top_bound = math.ceil(1.005 * x_start)  # bounds of sieve
    top_bound = int(math.ceil(n**0.5))+2*B

    # p_set = sieve_era(B)
    orig = [i*i-n for i in range(x_start, top_bound)]  # preserved list
    l = [math.log(i*i-n) for i in range(x_start,top_bound)]  # list of all the corresponding x^2-n that will be modified

#    b_smooth_list = []
    for p in p_set:
        #print(p)
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
                    print(str(x) + " " + str(orig[x-x_start]))
                    yield (x, orig[x-x_start])
                    # print(orig[j-x_start])
                #if len(b_smooth_list)>pi(B):
                    # return b_smooth_list

        elif len(sols) == 1:  # p | n
            print(str(p) + " " + str(orig[x-x_start]))
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
                    print(str(x) + " " + str(orig[x-x_start]))
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


def roots_congruent_mod_n(n, square, p_set, exponent_vector):
    pass


def gcd(a, b):
    return gcd(b, a % b) if a % b else b


def miller_rabin(n):
    m = n-1
    k = 0
    while m % 2 == 0:
        k += 1
        m = m / 2
    a = 2
    a = (a ** m) % n
    if a == 1 or a == n-1:
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


def factor(n, x, y):
    # TODO - return gcd(x-y, n) and n/gcd(x-y,n)
    return [0, 0]


def find_bound(n):
    return int(math.exp((2**(-0.5)*(math.log(2*n)*math.log(math.log(2*n)))**0.5)))


# testing with n
#n = 1127239
#B = int(math.exp((2**(-0.5)*(math.log(2*n)*math.log(math.log(2*n)))**0.5)))
#print(B)
#print(pi(B))

#sieve_quad_poly(n, B)
