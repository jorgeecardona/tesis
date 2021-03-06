#!/usr/bin/python

from sympy import *

def Symbolic_Matrix(name,size):
    M = zeros(size)
    for i in xrange(M.lines):
        for j in xrange(M.cols):
            M[i,j] = Symbol("%s[%i,%i]"%(name,i,j))
    return M


def luenberger(A, C, Polos):

    
    if reduce(lambda a,b: a and (type(b) is list),Polos,True):
        return map(lambda i: luenberger(A,C,i), Polos)

    if A.lines != A.cols:
        raise Exception, "size A mismatch"

    if A.lines != C.cols:
        raise Exception, "size C mismatch"

    K = Symbolic_Matrix("K", (A.lines, 1))

    x = var("x")
    p1 = (A - K*C).berkowitz_charpoly(x)

    p2 = reduce(lambda p, polo: p*(x-polo), Polos, Poly(1,x))
    
    eqs = map(lambda i: p1.coeff(i) - p2.coeff(i), xrange(p1.degree))
    
    sol = solve(eqs, K.mat)

    K.mat = map(lambda i: sol[i], K.mat)
    
    return K
    
