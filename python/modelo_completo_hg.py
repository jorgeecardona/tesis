from hg import hg
from sympy import *


# Create a linear object
sys = hg()

# Define global parameters
m_1, m, Y, A, l_1, l_2, l_3, l_4, l_5, r, k_f, L_k, rho_o, rho_e, g, I, l_1_, theta_ = sys.parameters("m_1 m Y A l_1 l_2 l_3 l_4 l_5 r k_f L_k rho_o rho_e 9.8 I l_1_ theta_")

# Define the sizes of the state and the input
x = sys.state(3)
u = sys.input(2)

# Define global functions that can be defined in an explicit way, i.e.
# y = f(x_1,x_2,...x_n) \forall n \in N and x_i \notequal y
def L(beta_1, beta_2, theta):
    return L_k + 2 * r * (beta_1 + 1/tan(beta_1) + beta_2 + 1/tan(beta_2)) + (1/sin(beta_1) + 1/sin(beta_2)) * (l_4 - l_1_ * sin(theta + theta_))

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
    lambda y,x: (l_3 - l_2) * sin(y) - l_4 * cos(y), # f(y, x_1, x_2, x_3. ... , x_n)
    lambda y,x: -l_1_ * sin(y + x + theta_) + 2 * r, # g(y, x_1, x_2, x_3, ... , x_n)
    [(x[0],)]
    )

beta_2 = sys.function(
    "beta_2",
    lambda y,x: l_2 * sin(y) - l_4 * cos(y), # f(y, x_1, x_2, x_3. ... , x_n)
    lambda y,x: l_1_ * sin(y - x - theta_) + 2 * r, # g(y, x_1, x_2, x_3, ... , x_n)
    [(x[0],)]
    )

f = lambda theta, M: T(rho(L(beta_1(theta),beta_2(theta),theta),M)) * (sin(theta + theta_ + beta_1(theta)) - sin(theta + theta_ - beta_2(theta)))

h = lambda theta, M: rho(L(beta_1(theta), beta_2(theta), theta), M)

# State function
sys.f(Matrix(
        [x[1],
         (l_1_ * f(x[0],x[2]) - l_1_ * m * g * cos(x[0] + theta_) - k_f * x[1]) / (I + l_1_**2 * m_1),
         0
         ])
      )

sys.g(Matrix(
        [[0,0],
         [0,0],
         [A * r * rho_e , - A * r * h(x[0], x[2])]]))

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

print "Para modelica:"

for i in xrange(J_inv.lines):
    for j in xrange(J_inv.cols):
        print "    Jac_inv[%d,%d] = "%(i+1,j+1),str(J_inv[i,j]).replace("beta_1(X_e[1,1])", "core.beta_1_aux").replace("beta_2(X_e[1,1])", "beta_2").replace("**", "^"), ";"

for i in xrange(J.lines):
    for j in xrange(J.cols):
        print "    Jac[%d,%d] = "%(i+1,j+1),str(J[i,j]).replace("beta_1(X_e[1,1])", "beta_1").replace("beta_2(X_e[1,1])", "beta_2").replace("**", "^"), ";"


    
from sympy import I
print sym.compute_h([-10 + 4 * I,-10 - 4*I,-10])
