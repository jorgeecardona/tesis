
def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 1233
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        
        else:
            return
        
        indices[i] += 1
        
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        
        yield tuple(pool[i] for i in indices)


def split(P):
    pass

def minimum(P):

    n = len(P)

    if n == 0:
        raise Exception, "No points."

    if n == 1:
        return P[0]

    m = len(P[0])

    if any([len(p) != m for p in P]):
        raise Exception, "Inconsistent point sizes."

    if m == 1:
        if any([(a * b < 0) for a, b in combinations(P,2)]):
            return array([0])
        
        return P[array([abs(p) for p in P]).argmin()]

    if m == n:
        # base solver
        pass
    
    R = split(P)

    if (len(R[1]) == 0) :
        #solve ()

        pass
    else:
        # Minimize each polytope
        R = [minimum(r) for r in R]
        
        # Select the minimum
        return R[array([dot(r,r) for r in R]).argmin()]
