from linear import linear
from sympy import *
from L import calcular_L


# Define global parameters
m_1     = 0.860  # Masa de la polea
m_2     = 0.407  # Masa del pivote
Y       = 2e6    # Modulo de Young de la cinta
A       = 14e-6  # Area transversal de la cinta
l_1     = 0.29   # Ver Diagramas del modelo
l_2     = 0.38   # Ver Diagramas del modelo
l_3     = 0.41   # Ver Diagramas del modelo
l_4     = 0.50   # Ver Diagramas del modelo
r       = 0.05   # Radio de las poleas
k_f     = 0.1    # Constante de friccion
L_k     = 0.3    # Ver Diagramas del modelo
rho_o   = 714    # Densidad de la cinta en reposo
rho_e   = 714    # Densidad de la cinta por fuera del sistema
g       = 9.81   # Gravedad
I       = 0.012  # Inercia del pivote

# Create a linear object, and define sizes
sys = linear()

x = sys.state(3)
u = sys.input(2)

# Define global functions that can be defined in an explicit way, i.e.
# y = f(x_1,x_2,...x_n) \forall n \in N and x_i \notequal y
def L(beta_1, beta_2, theta):
    return L_k + (
        2 * r * (beta_1 + 1/tan(beta_1) + beta_2 + 1/tan(beta_2))) + (
        (1/sin(beta_1) + 1/sin(beta_2)) * (l_4 - l_1 * sin(theta)))

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
    lambda y,x: -l_1 * sin(y + x) + 2 * r, # g*(y,x_1,x_2,... ,x_n)
    [(x[0],)]
    )

beta_2 = sys.function(
    "beta_2",
    lambda y,x: l_2 * sin(y) - l_4 * cos(y), # f*(y,x_1,x_2,...,x_n)
    lambda y,x: l_1 * sin(y - x) + 2 * r, # g*(y,x_1,x_2,...,x_n)
    [(x[0],)]
    )

f = lambda theta, M: T(rho(L(beta_1(theta),beta_2(theta),theta),M)) * (
    sin(theta + beta_1(theta)) - sin(theta - beta_2(theta)))

h = lambda theta, M: rho(L(beta_1(theta), beta_2(theta), theta), M)

# State function and output function
sys.f(
    Matrix(
        [x[1],

         (l_1 * f(x[0],x[2]) - l_1 * (m_1 + m_2) * g * cos(x[0]) - (
                    k_f * x[1])) / (I + l_1 * m_1 * r),

         A * r * (rho_e * u[0] - h(x[0],x[2]) * u[1])
         ])
    )

sys.h(Matrix([x[0]]))

# Linearization!!!
A,B,C = sys.linearize(Matrix([0, 0, 0.01389]), Matrix([5, 5.9699]))

print "A=\n", A
print "B=\n", B
print "C=\n", C

# print latex(A)
# print latex(B)
# print latex(C)


Raices = []
Raices.append([-0.5 + 1.0*I, -0.5 - 1.0*I, -2.5])
Raices.append([-2.5 + 5.0*I, -2.5 - 5.0*I, -12.5])
Raices.append([-12.5 + 25.0*I, -12.5 - 25.0*I, -62.5])

# for Raiz in Raices:

#     L =  calcular_L(A,C,Raiz)

#     print "*********************   R   **********************"
#     print Raiz

#     print "*********************   L   **********************"
#     print latex(L.evalf())

#     print "*********************   A   **********************"
#     print A.evalf()

#     print "*********************   B   **********************"
#     print B.row_join(L).evalf()

#     print "*********************   C   **********************"
#     print C.evalf()

#     print "******************** SCICOS  *********************"

from luenberger import luenberger
print luenberger(A, C, Raices)
