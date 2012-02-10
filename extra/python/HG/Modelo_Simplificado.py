from HG import HG, calcular_H
from sympy import *


# Create a linear object
sys = HG()

# Define global parameters
M_p, Y, A, l, r, k_f, L_k, rho_o, rho_e, g = sys.parameters("M_p Y A l r k_f L_k rho_o rho_e g")

# Define global functions that can be defined in an explicit way, i.e.
# y = f(x_1,x_2,...x_n) \forall n \in N and x_i \notequal y
def L(beta):
    return L_k + 4 * r * ( beta - tan(beta)) + l / cos(beta)


def T(rho):
    return Y * A * (rho_o /rho - 1)

def rho(L,M):
    return M /(A * L)

L = Function("L")
T = Function("T")
rho = Function("rho")

# Define the sizes of the state and the input
x = sys.state(3)
u = sys.input(2)

# Define functions that can't be defined as functions above, i.e.
# f(y,x_1,x_2,x_3,...,x_n) = g(y,x_1,x_2,x_3,...,x_n)
# The method receive the name of the function, the f(y,x_1,x_2,...,x_n) 
# function, the g(y,x_1,x_2,...,x_n) function, and a list with tuples 
# representing the args that call the function (this is a hack, 
# because the actual matching system is not as good as i want).

beta = sys.function(
    "beta",
    lambda y,x: l * sin(y) - 4*r, # f(y, x_1, x_2, x_3. ... , x_n)
    lambda y,x: 2 * x * cos(y), # g(y, x_1, x_2, x_3, ... , x_n)
    [(x[0],)]
    )

# State function
sys.f(Matrix(
        [x[1],
         g - 2 * T(rho(L(beta(x[0])),x[2])) * sin(beta(x[0]))/M_p - k_f * x[1] / M_p,
         A * r * (rho_e * u[0] - rho(L(beta(x[0])),x[2]) * u[1] )
         ])
      )

# Output function
sys.h(Matrix([x[0]]))

# HG!!!
M, J, J_inv = sys.hg()

print "Mapa="
print(M)
print "\n\n"

print "Jacobiano="
print(J)
print "\n\n"

print "Inversa del Jacobiano="
print(J_inv)
print "\n\n"


I   = 0.023
pprint(calcular_H([-0.2+0.4*I,-0.2-0.4*I,-1]))
pprint(calcular_H([-2+4*I,-2-4*I,-10]))
pprint(calcular_H([-20+40*I,-20-40*I,-100]))
pprint(calcular_H([-200+400*I,-200-400*I,-1000]))
