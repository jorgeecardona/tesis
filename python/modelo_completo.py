from linear import linear, Symbolic_Matrix
from sympy import *
from luenberger import luenberger
from kalman import kalman

# Define global parameters
m_1     = 0.5     # Masa de la polea
m       = 1.013   # MAsa conjunta de la polea y el pivote.
Y       = 3.744e6 # Modulo de Young de la cinta
A       = 31e-6   # Area transversal de la cinta
l_1     = 0.27    # Ver Diagramas del modelo
l_2     = 0.384   # Ver Diagramas del modelo
l_3     = 0.41    # Ver Diagramas del modelo
l_4     = 0.50    # Ver Diagramas del modelo
l_5     = 0.0952
l_1_    = sqrt(l_1**2 + l_5**2)
theta_  = atan(l_5 / l_1)
r       = 0.05    # Radio de las poleas
k_f     = 0.89182 # Constante de friccion
L_k     = 2.5     # Ver Diagramas del modelo
rho_o   = 1215.7  # Densidad de la cinta en reposo
rho_e   = 1041    # Densidad de la cinta por fuera del sistema
g       = 9.81    # Gravedad
I       = 0.13972 # Inercia del pivote

# Create a linear object, and define sizes
sys = linear()

x = sys.state(3)
u = sys.input(2)

# Define global functions that can be defined in an explicit way, i.e.
# y = f(x_1,x_2,...x_n) \forall n \in N and x_i \notequal y
def L(beta_1, beta_2, theta):
    return L_k + (
        2 * r * (beta_1 + 1/tan(beta_1) + beta_2 + 1/tan(beta_2))) + (
        (1/sin(beta_1) + 1/sin(beta_2)) * (l_4 - l_1_ * sin(theta + theta_)))

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

beta_1 = sys.function(
    "beta_1",
    lambda y,x: (l_3 - l_2) * sin(y) - l_4 * cos(y), # f*(y,x_1,x_2,...,x_n)
    lambda y,x: -l_1_ * sin(y + x + theta_) + 2 * r, # g*(y,x_1,x_2,... ,x_n)
    [(x[0],)]
    )

beta_2 = sys.function(
    "beta_2",
    lambda y,x: l_2 * sin(y) - l_4 * cos(y), # f*(y,x_1,x_2,...,x_n)
    lambda y,x: l_1_ * sin(y - x - theta_) + 2 * r, # g*(y,x_1,x_2,...,x_n)
    [(x[0],)]
    )

f = lambda theta, M: T(rho(L(beta_1(theta),beta_2(theta),theta),M)) * (
    sin(theta + theta_ + beta_1(theta)) - sin(theta + theta_- beta_2(theta)))

h = lambda theta, M: rho(L(beta_1(theta), beta_2(theta), theta), M)

# State function and output function
sys.f(
    Matrix(
        [x[1],

         (l_1_ * f(x[0],x[2]) - l_1_ * m * g * cos(x[0] + theta_) - (
                    k_f * x[1])) / (I + l_1_**2 * m_1),

         A * r * (rho_e * u[0] - h(x[0],x[2]) * u[1])
         ])
    )

sys.h(Matrix([x[0]]))

# Linearization!!!
A,B,C = sys.linearize(Matrix([0.02824, 0, 0.13233]), Matrix([5.0, 4.46]))

print "A=\n", A
print "B=\n", B
print "C=\n", C

from sympy import I
K = luenberger(A, C, [-10 + 4 * I, -10 - 4 * I, -10])

print "K = \n", K

V_2 = Matrix([1e-2])
V_1 = kalman(A, C, K, V_2)

print "V_1 = \n", V_1
print "V_2 = \n", V_2
