from sympy import *
from linear import linear
from luenberger import luenberger


# Define global parameters
M_p     = 0.86   # Masa de la polea
Y       = 2e6    # Modulo de Young de la cinta
A       = 14e-6  # Area transversal de la cinta
l       = 0.41   # Distancia entre poleas 2 y 3
r       = 0.05   # Radio de las poleas
k_f     = 1.0    # Constante de friccion
L_k     = 0.3    # Ver Diagrama del modelo
rho_o   = 714    # Densidad de la cinta en reposo
rho_e   = 714    # Densidad de la cinta por fuera del sistema

# Create a linear object, and define sizes
sys = linear()

x = sys.state(3)
u = sys.input(2)

# Define global functions that can be defined in an explicit way, i.e.
# y = f(x_1,x_2,...x_n) \forall n \in N and x_i \notequal y
def L(beta):
    return L_k + 4 * r * ( beta - tan(beta)) + l / cos(beta)

def T(rho):
    return Y * A * (rho_o /rho - 1)

def rho(L,M):
    return M /(A * L)

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

f = lambda b, M: T(rho(L(beta(b)),M)) * sin(beta(b))

h = lambda b, M: rho(L(beta(b)),M)

# State function and output function
sys.f(Matrix(
        [x[1],

         9.81 - 2 * f(x[0], x[2])/M_p - k_f * x[1] / M_p,

         A * r * (rho_e * u[0] - h(x[0], x[2]) * u[1] )
         ])
      )

sys.h(Matrix([x[0]]))

# Linearization!!!
A,B,C = sys.linearize(Matrix([0.5,0,0.01417]), Matrix([5, 5.7691]))

print "A=\n", A
print "B=\n", B
print "C=\n", C

K1, K2, K3 = luenberger(A,C,[
        [-0.1 + 0.2*I, -0.1 - 0.2*I, -0.5],
        [-0.5 + 1.0*I, -0.5 - 1.0*I, -2.5],
        [-2.5 + 5.0*I, -2.5 - 5.0*I, -12.5]])

print "K1=\n",K1
print "K2=\n",K2
print "K3=\n",K3
