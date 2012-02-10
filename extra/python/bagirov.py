from numpy import array, matrix, sin, sqrt, dot, cos, ix_, zeros, concatenate, abs, log10, exp, ones
from numpy.random import random
from numpy.linalg import norm
from itertools import takewhile, count

def discret_gradient(f, x, g, i=None, E = None, z=sin, lam=0.1, alpha=0.1):

    # Space size
    n = len(x)

    if E is None:
        E = ones(n)

    Ej = [array([(e * alpha**(k+1) if k<=j else 0) for k,e in zip(range(n),E)]) for j in xrange(n)]

    # Calculate points
    X = [x + lam * g] + [x + lam * g + z(lam) * e for e in Ej] + [x]

    # Evaluate points
    F = [f(x) for x in X]

    # Get argmax
    if i is None:
        i = abs(g).argmax()

    # Calculate discret gradient
    gradient = array([float(F[j+1] - F[j]) / (e * z(lam) * alpha**(j+1)) if (j!=i) else 0 for j,e in zip(xrange(n),E)])

    gradient[i] = float(F[0] - F[-1] - lam * (gradient * g).sum()) / (lam*g[i])

    return gradient

def find_min_point(P):
    print "Calling find_min with P: ", P

    if len(P) == 1:
        return P[0]

#    print "Called with: ", P[:10]

    # Step 0. Choose a point from C(P)
    x  = P[array([dot(p,p) for p in P]).argmin()]

    while True:

        # Step 1. \alpha_k := min{x_{k-1}^T p | p \in P}        
        Pk = [P[array([dot(x,p) for p in P]).argmin()]]

        alpha = dot(x,Pk[0])

        if (dot(x,x) <= alpha):
            return x 

        # Step 2. P_k := { p | p \in P and x_{k-1}^T p = \alpha_k}    
        Pk = Pk + [p for p in P if dot(x,p) == alpha and not (p == Pk[0]).all()]

        # Check if P == P_
        if array([array([(p == q).all() for q in Pk]).any() for p in P]).all():
            return x

        y = find_min_point(Pk)
        
        # Step 3. B_k := min{y_k^T p | p  \in P\P_k}
        beta = array([dot(y,p) for p in P if not array([(p == q).all() for q in Pk]).any()]).min()
        
        if (dot(y,y) <= beta):
            return y

        # Step 4. 
#        print dot(x,p-y)
        lam = array([float(dot(x,p-y)) / dot(y-x,y-p) for p in P if not array([(p==q).all() for q in Pk]).any() and dot(y-x,y-p) > 0]).min()
        

#        print lam
        if lam < 1e-15:
            return x

        x += lam * (y-x)


    raise Exception, "Error"



def descent_direction(f, x, E=None, c=0.1, lam=0.1, delta=0.1, alpha=0.1, f_ = None):

    # Space size
    n = len(x)

    if E is None:
        E = ones(n)

    # initial direction
    g = 2*random(n)-1
    g = g / norm(g)

    D = []
    while True:
        
        v = discret_gradient(
            f,
            x, 
            g, 
            lam=lam, 
            alpha=alpha, 
            E = E,
            )


        if (len(D)>0) and (v == D[-1]).all():
            return g, W

        print "Gradient: ", v
        D.append(v)
        print "D: ", D

#        print D  , g
        W = find_min_point(D)
        print "W: ", W

        if norm(W) <=  delta:
            return g, W

        g = -W/norm(W)

        if f(x + lam * g) - f(x) <= - c * lam * norm(W):
            return g, W

    raise Exception, "Error"

def discret_gradient_method(f, xo, max_calls = 100000, E = None, c1=1e-2, c2=1e-2, lam=1e-2, delta=1e-4, beta = 1, alpha=1e-2, args=[], kwords={}, use_map = False):

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
                raise Exception, "Maximum number of function call."
            
            if self.use_map:
                if tuple(X) in self.map:
                    return  self.map[tuple(X)]

            res = self.f(X,*self.args,**self.kwords)

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

    c2 = c1 * c2
    try:
        while True:
            lam  = lam_0 * beta ** len(X)
            
            f_ = fun(X[-1][-1])
            
            g,v = descent_direction(
                fun, 
                X[-1][-1], 
                E = E,
                c = c1, 
                lam = lam,
                delta = delta, 
                alpha = alpha, 
                f_ = f_,
                )
            
            if norm(v) <= delta:
                X.append([X[-1][-1]])
                continue
            
            v_ = norm(v)
            
            sigma = list(
                takewhile(
                    lambda t: fun(X[-1][-1] + t*lam*g) - f_ <= - c2*t*lam*v_, 
                    (exp(i) for i in count())
                    )
                )
            
            if len(sigma) == 0:
                return X
            
            sigma = sigma[-1]*lam
            
            X[-1].append(X[-1][-1] + sigma * g)

            print "X: ", X[-1][-1]

    except:
        return X

            
            


    #Step 1: Choose any g^1 \in S_1

#print "Min: ", find_min_point([array([-2,1]), array([3,0]),array([0,2])])

c = 0

def f(x):
    global c 
    c +=1
#    print "Eval %d :"%(c), x
    res = cos(x[0]) + sin(float(x[1])/1000)
#    print "Eval: ", res
    return res


print " , ".join(["%0.12f"%(i) for i in discret_gradient_method(
    f,
    array([2,200]),
    E = array([1,1000]),
    delta = 1e-8,
    c1 = 1e-2,
    lam = 1e-2,
    beta = (0.01)**(1/float(1000)),
    max_calls = 10000,
    )[-1][-1]])

print "Eval: ", c
