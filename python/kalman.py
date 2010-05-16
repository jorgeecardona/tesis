from linear import Symbolic_Matrix
from sympy import *

def kalman(A, C, K, V_2):
    
    P_sym = Symbolic_Matrix("p",(3,3))

    # Usar ecuacion de Riccati para encontrar el valor de P
    res = solvers.solve(tuple((P_sym * C.transpose() * V_2.inv() - K).mat), P_sym.mat)
    P = zeros(P_sym.shape)
    P.mat = [p if not p in res else res[p] for p in P_sym.mat]

    # Usar ecuacion de Riccati para encontrar el valor de V_1
    V_1 = - (P * A.transpose() + A * P - P * C.transpose() * V_2.inv() * C * P)    
    res = solvers.solve(tuple([V_1[i,j] - V_1[j,i] for i in xrange(3) for j in xrange(3) if i < j] + [V_1[0,1], V_1[0,2], V_1[1,2]]), P_sym.mat)
    P.mat = [p if not p in res else res[p] for p in P.mat]
    V_1 = - (P * A.transpose() + A * P - P * C.transpose() * V_2.inv() * C * P)
    
    return V_1
