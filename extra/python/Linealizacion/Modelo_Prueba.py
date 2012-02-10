from Linear import Linear
from L import calcular_L
from sympy import *


# Define global parameters
A_1, A_c, A_2, K_1, K_c, K_2, b_c = var("A_1 A_c A_2 K_1 K_c K_2 b_c")

# Create a linear object
sys = Linear()

# Define the sizes of the state and the input
x = sys.state(3)
u = sys.input(1)

# Define functions that can't be defined as functions above, i.e.
# f(y,x_1,x_2,x_3,...,x_n) = g(y,x_1,x_2,x_3,...,x_n)
# The method receive the name of the function, the f(y,x_1,x_2,...,x_n) 
# function, the g(y,x_1,x_2,...,x_n) function, and a list with tuples 
# representing the args that call the function (this is a hack, 
# because the actual matching system is not as good as i want).

# State function
sys.f(Matrix(
        [u[0]/A_1 -K_1 / A_1 * pow(x[0], 2.475),
         (K_1 / A_c) * pow(x[0], 2.475) - (K_c * b_c / A_c) * pow(x[1], 1.8),
         (K_c * b_c / A_2) * pow(x[1], 1.8) - (K_2 / A_2) * x[2]
         ])
      )

# Output function
sys.h(Matrix([x[2]]))

# Linearization!!!
A,B,C = sys.linearize(x, u)


print "A ="
print(A)

print "B ="
print(B)

print "C ="
print(C)


