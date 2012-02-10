from linear import linear
from L import calcular_L
from sympy import *


# Define global parameters
M_p     = 0.86
Y       = 2e6
A       = 14e-6
l       = 0.41
r       = 0.05
k_f     = 0.1
L_k     = 0.3
rho_o   = 714
rho_e   = 714

# Define global functions that can be defined in an explicit way, i.e.
# y = f(x_1,x_2,...x_n) \forall n \in N and x_i \notequal y
def L(beta):
    return L_k + 4 * r * ( beta - tan(beta)) + l / cos(beta)

def T(rho):
    return Y * A * (rho_o /rho - 1)

def rho(L,M):
    return M /(A * L)


# Create a linear object
sys = linear()

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
         9.81 - 2 * T(rho(L(beta(x[0])),x[2])) * sin(beta(x[0]))/M_p - k_f * x[1] / M_p,
         A * r * (rho_e * u[0] - rho(L(beta(x[0])),x[2]) * u[1] )
         ])
      )

# Output function
sys.h(Matrix([x[0]]))

# Linearization!!!
A,B,C = sys.linearize(Matrix([0.5,0,0.01417]), Matrix([5, 5.7690]))


print "A ="
print(A)

print "B ="
print(B)

print "C ="
print(C)

print latex(A)
print latex(B)
print latex(C)

Raices = []
Raices.append([-0.1 + 0.2*I, -0.1 - 0.2*I, -0.5])
Raices.append([-0.5 + 1.0*I, -0.5 - 1.0*I, -2.5])
Raices.append([-2.5 + 5.0*I, -2.5 - 5.0*I, -12.5])

for Raiz in Raices:

    L =  calcular_L(A,C,Raiz)

    print "*********************   R   **********************"
    print Raiz

    print "*********************   L   **********************"
    print latex(L.evalf())

    print "*********************   A   **********************"
    print A.evalf()

    print "*********************   B   **********************"
    print B.row_join(L).evalf()

    print "*********************   C   **********************"
    print C.evalf()

    print "******************** SCICOS  *********************"

