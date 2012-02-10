#!/usr/bin/python

from sympy import *
from random import random

def print_scicos(M):
    print "[",
    for i in xrange(M.lines):
        for j in xrange(M.cols):
            print M[i,j],
            if j < M.cols-1:
                print " ",
        if i < M.lines - 1:
            print ";",
    print "]"


def calcular_L(Matriz_A,Matriz_C,Raices):


    P = 0
    for i in Raices:
        P += abs(i)/len(Raices)


    # Definir simbolicamente los elementos de L
    L_elementos = []
    L_elementos_raw = []
    for i in xrange(Matriz_A.cols):
        L_elementos.append([])
        for j in xrange(Matriz_C.lines):
            if j > 0:
                L_elementos[-1].append(P*random())
            else:
                L_elementos[-1].append(var("L_%d"%(i)))
                L_elementos_raw.append(L_elementos[-1][-1])
    
    # Definir L
    L = Matrix(L_elementos)
    
    # Crear el polinomio a partir del determninante.
    var('lam')
    p = Poly((lam*eye(Matriz_A.cols)-(Matriz_A-L*Matriz_C)).berkowitz_det(),lam)

    # Crear el otro polinomio a partir de las raices.
    p1 = 1
    for i in Raices:
        p1 *=Poly(lam-i,lam)
    
    # Comparar grados
    if p.degree != p1.degree:
        raise Exception, "Polynomial degree mismatch"
    
    # Crear sistema de ecuaciones
    sys = ()
    for i in xrange(p.degree + 1):
        sys += ((p.coeff(i)-p1.coeff(i)),)
   
    r = solve(sys,L_elementos_raw,simplified=False)

    # Reconstruir L con soluciones
    for i in xrange(L.lines):
        L[i,0] = r[L[i,0]]

    # Retornar L
    return L
