from numpy import array, matrix, sin, sqrt, dot, cos, ix_, zeros, concatenate, abs, log10, exp, ones
from numpy.random import random
from numpy.linalg import norm
from mpmath import mpf, mp, almosteq
mp.dps=100
from itertools import takewhile, count

from wolfe import minimum

def discret_gradient(f, x, g, i=None, z=sin, lam=0.1, alpha=0.1):

    # Space size
    n = len(x)
    e_j =lambda alpha, j: alpha**(1.0 + j * 0.1)


    E = [array([(e_j(alpha,k+1) if k<=j else 0.0) for k in range(n)]) for j in xrange(n)]

    # Calculate points
    X = [x + lam * g] + [x + lam * g + z(lam) * e for e in E] + [x]
    print "X: ", X

    # Evaluate points
    F = [f(x) for x in X]
    print "F: ", F

    # Get argmax
    if i is None:
        i = abs(g).argmax()

    # Calculate discret gradient
    gradient = array([float(F[j+1] - F[j]) / (z(lam) * e_j(alpha,j+1)) if (j!=i) else 0.0 for j in xrange(n)])
    gradient[i] = float(F[0] - F[-1] - lam * (gradient * g).sum()) / ( lam * g[i] )

    return gradient

def descent_direction(f, x, c=0.1, lam=0.1, delta=0.1, alpha=0.1 ):

    # Space size
    n = len(x)

    # initial direction
    g = 2.0 * random(n) - 1.0
    g = g / norm(g)

    D = []
    while True:
        
        v = discret_gradient(
            f,
            x, 
            g, 
            lam=lam, 
            alpha=alpha, 
            )

        

        print "Gradient: ", v
        D.append(v)
        print "D: ", D

        W = minimum(D)

        print "W: ", W
        print "|W|: ", norm(W)

        if (len(D)>1) and (norm(W) >= norm(minimum(D[:-1]))):
            print "Error"

            W = (delta**2) * g / norm(g)

            return -W/norm(W) , W

        if norm(W) <=  delta:
            return g, W

        g = -W/norm(W)
        
        if f(x + lam * g) - f(x) <= - c * lam * norm(W):
            return g, W

        print "c < ", (f(x + lam * g) - f(x)) / ( - lam * norm(W))

    raise Exception, "Error"

class MaxCalls(Exception):
    pass

def discret_gradient_method(f, xo, max_calls = 100000, E = None, c1=1e-2, c2=1e-2, lam=1e-2, delta=1e-4, beta = 1, alpha=1e-2, args=(), kwords={}, use_map = False):

    class Fun:
        def __init__(self, f, max_calls = 10000, use_map = False, args = (), kwords = {}):
            self.args = args
            self.kwords = kwords
            self.f = f
            self.counter = 0
            self.max_calls = max_calls
            self.use_map = use_map
            self.map ={}
            
        def __call__(self, X):

            self.counter += 1
            if self.counter >= self.max_calls:
                raise MaxCalls
            
            if self.use_map:
                if tuple(X) in self.map:
                    return  self.map[tuple(X)]

            res = self.f(X, *self.args, **self.kwords)

            if self.use_map:
                self.map[tuple(X)] = res

            return res

    fun = Fun(
        f, 
        max_calls = max_calls,
        args = args,
        kwords = kwords,
        use_map = use_map
        )

    X = [[array([float(x) for x in xo])]]

    lam_0 = lam
    alpha_0 = alpha

    c2 = c1 * c2
    try:
        while True:
            lam  = lam_0 * beta ** (len(X)-1)
            alpha = alpha_0 * beta ** (len(X)-1)
            
            print "Step: ", len(X)
            print " lambda: ", lam
            print " alpha: ", alpha

            g,v = descent_direction(
                fun, 
                X[-1][-1], 
                c = c1, 
                lam = lam,
                delta = delta, 
                alpha = alpha, 
                )

            print "Direction: ", g

            if norm(v)  == 0:
                return X[-1][-1], X, fun.counter

            if norm(v) <= delta:
                X.append([X[-1][-1]])
                continue

            sigma = list(
                takewhile(
                    lambda t: fun(X[-1][-1] + t*lam*g) - fun(X[-1][-1]) <= - c2 * t * lam * norm(v), 
                    (i for i in count())
                    )
                )
            
            if len(sigma) == 0:
                print "f(X): ", fun(X[-1][-1])
                print "f(C + lam*g): ", fun(X[-1][-1] + lam * g)
                print "c1 "
                return X[-1][-1], X, fun.counter
            
            sigma = sigma[-1]*lam
            
            print "Sigma: ", sigma

            X[-1].append(X[-1][-1] + sigma * g)

            print "X: ", X[-1][-1]

    except MaxCalls:
        print "Maximum number of function calls."
        return X[-1][-1], X, fun.counter

            
            


    #Step 1: Choose any g^1 \in S_1

#print "Min: ", find_min_point([array([-2,1]), array([3,0]),array([0,2])])

if __name__ == "__main__":
    def f(x):
        return sin(float(x[0])/1000)

    res = discret_gradient_method(
        lambda x: (1-x[0])**2 + 100*(x[1] - x[0]**2)**2,
        array([-1,1]),
        delta = 1e-8,
        c1 = 1e-2,
        c2 = 1e-2,
        lam = 1e-1,
        alpha = 1e-2,
        beta = 0.8,
        max_calls = 8000,
        )
    
    print "X = ", res[0]
    print "Evaluations: ", res[2]
    print "Iterations: ", len(res[1])
    print "Iterations: ", sum(len(r) for r in res[1])
    
    
