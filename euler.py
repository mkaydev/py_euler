import math
import queue

def concat_numbers(numbers_seq):
    concat_str = "".join(str(n) for n in numbers_seq)
    return int(concat_str)

def get_digit_sum(n):
    return sum(int(c) for c in str(n))

def is_pos_int(n, error=0.000000001):    
    return n > 0 and abs(round(n) - n) <= error

def is_tri_num(tn):
    n1 = -1/2 + math.sqrt(1/4 + 2 * tn)
    n2 = -1/2 - math.sqrt(1/4 + 2 * tn)
    n = max(n1, n2)
    return is_pos_int(n)

def is_squ_num(sn):
    if sn < 1: return False
    return is_pos_int(math.sqrt(sn))

def is_pent_num(pn):
    n1 = 1/6 + math.sqrt(1/36 + 2/3 * pn)
    n2 = 1/6 - math.sqrt(1/36 + 2/3 * pn)
    n = max(n1, n2)
    return is_pos_int(n)

def is_hex_num(hn):
    n1 = 1/4 + math.sqrt(1/16 + 1/2 * hn)
    n2 = 1/4 - math.sqrt(1/16 + 1/2 * hn)
    n = max(n1, n2)
    return is_pos_int(n)

def is_hept_num(hn):
    n1 = 3/10 + math.sqrt(9/100 + 2/5 * hn)
    n2 = 3/10 - math.sqrt(9/100 + 2/5 * hn)
    n = max(n1, n2)
    return is_pos_int(n)

def is_oct_num(on):
    n1 = 1/3 + math.sqrt(1/9 + 1/3 * on)
    n2 = 1/3 - math.sqrt(1/9 + 1/3 * on)
    n = max(n1, n2)
    return is_pos_int(n)

def get_digit_count(n, base=10):
    assert n > 0 and base > 1
    return math.floor(math.log(n, base)) + 1

# Euclidean algorithm
# https://en.wikipedia.org/wiki/Euclidean_algorithm
def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

# https://en.wikipedia.org/wiki/Least_common_multiple#Reduction_by_the_greatest_common_divisor
def lcm(a, b):
    return (abs(a) / gcd(a, b)) * abs(b) 

def nCr(n,r):
    f = math.factorial
    return (f(n) / f(r)) / f(n-r)

def is_palin(n):
    n_str = str(n)
    left = 0
    right = len(n_str) - 1

    while left < right:
        if n_str[left] != n_str[right]: return False

        left += 1
        right -= 1 
    return True

# heron's formula
# https://en.wikipedia.org/wiki/Heron%27s_formula
def get_triangle_area(a, b, c):
    assert a > 0 and b > 0 and c > 0
    s = (a + b + c) / 2
    A = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return A  

def is_prime(n):
    if n == 1: return False
    if n == 2: return True
    if n % 2 == 0: return False

    for i in range(3, math.ceil(math.sqrt(n)), 2):
        if (n % i == 0): return False
    return True

# sieve of eratosthenes
# https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
def get_primes(n):
    is_prime = {i: True for i in range(2, n)}

    p = 2
    while p * p < n:
        for i in range(p * p, n, p):
            is_prime[i] = False

        found_next = False
        for i in range(p + 1, n):
            if is_prime[i]:
                p = i
                found_next = True
                break

        if not found_next: break

    primes = [i for i in is_prime if is_prime[i]]
    primes.sort()
    return primes

# Fermat's factorization method - n must be odd
# https://en.wikipedia.org/wiki/Fermat%27s_factorization_method
def get_prime_factor(n):
    assert n % 2 == 1
    
    a = int(math.ceil(math.sqrt(n)))
    b2 = a * a - n

    while int(math.sqrt(b2))**2 != b2:
        a += 1
        b2 = a * a - n

    return int(a - math.sqrt(b2))

def get_prime_factors(n):
    factors = {}
    if n == 1: return {}

    while n % 2 == 0:
        n = n / 2
        if 2 in factors:
            factors[2] += 1
        else:
            factors[2] = 1
    
    f = get_prime_factor(n)
    new_n = n

    while f != 1:
        ff = get_prime_factor(f)

        while ff != 1:
            f = ff
            ff = get_prime_factor(f)
        
        if f in factors:
            factors[f] += 1
        else:
            factors[f] = 1
        
        new_n = new_n / f
        f = get_prime_factor(new_n)

    f = int(new_n)
    if f in factors:
        factors[f] += 1
    else:
        factors[f] = 1

    if 1 in factors: factors.pop(1)
    return factors

# min (primitive) positive solution to positive pell equation x * x - D * y * y = 1
# using continued fraction expansion http://mathworld.wolfram.com/PellEquation.html
def get_pos_pell_solution(D):
    assert D > 1
    
    m = q1 = 0
    d = q = p1 = 1
    a = a0 = p = math.floor(math.sqrt(D))
    
    while p * p - D * q * q != 1:
        m = int(d * a - m) 
        d = int((D - m * m) / d) 
        a = int(math.floor((a0 + m) / d))

        p2, p1 = p1, p
        q2, q1 = q1, q

        p = a * p1 + p2
        q = a * q1 + q2
    return p, q

class Graph:
    def __init__(self, nodes, roots=None, ends=None):
        self.nodes = nodes
        self.roots = roots
        self.ends = ends

class GraphNode:
    def __init__(self, value):
        self.value = value
        self.children = []

    def add_connection(self, child_node):
        if not child_node in self.children:
            self.children.append(child_node)
        
    def __lt__(self, other):
        return self.value < other.value

    def __gt__(self, other):
        return self.value > other.value

    def __le__(self, other):
        return self.value <= other.value

    def __ge__(self, other):
        return self.value >= other.value

# djikstra's algorithm for minimum path from source to a target within a graph
# https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
def get_djikstra_min_path(graph, source, targets):
    dist = {}
    visited = {}
    previous = {}
    vertices = graph.nodes

    for v in vertices:
        v.dist = None
        v.previous = None
        v.visited = False

    source.dist = source.value
    Q = queue.PriorityQueue()
    Q.put((source.dist, source))

    while not Q.empty():
        dist_u, u = Q.get()
        
        # note to self: in-operator uses ==-operator, NOT is-operator
        if u in targets: 
            target = u
            break
        u.visited = True

        neighbors = u.children

        for v in neighbors:
            alt = dist_u + v.value
            
            if v.dist is None or alt < v.dist:
                v.dist = alt
                v.previous = u

                if not v.visited: Q.put((alt, v))

    S = []
    u = target

    while u.previous is not None:
        S.append(u.value)
        u = u.previous

    assert u is source
    S.append(u.value)
    S.reverse()
    
    return S