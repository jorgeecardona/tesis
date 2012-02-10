#!/usr/bin/python


from HG import *
from sympy import Matrix, Function, cos, sin, tan, diff, pprint
import scipy, scipy.optimize

cot = lambda x: 1/tan(x)
csc = lambda x: 1/sin(x)

m1,m2,I, k_f, r, l1, l2, l3, l4, Y, A  = var("m1 m2 I k_f r l1 l2 l3 l4 Y A")

var("t")

x = Symbolic_Matrix("x",(3,1))
u = Symbolic_Matrix("u",(2,1))

X = Matrix([ -0.45230248276250362, 0 , 4.9132455972876734])
U = Matrix([ 0.3600036000363816, 0.3600036000363816])



f = zeros((3,1))
beta_1   = Function("beta1")(x[0])
beta_2   = Function("beta2")(x[0])
db1 = l1 * cos (beta_1)**2 * x[1] * (cos(x[0]) - tan(beta_1) * sin(x[0]) )/(2 * r * sin(beta_1) - l1 * cos(x[0]) + l2 - l3 )
db2 = l1 * cos (beta_2)**2 * x[1] * (cos(x[0]) + tan(beta_2) * sin(x[0]) )/(2 * r * sin(beta_2) + l1 * cos(x[0]) - l2 )
L   = 2 * r * (beta_1 + cot(beta_1) + beta_2 + cot(beta_2)) + (csc(beta_1) + csc(beta_2)) * (l4 - l1 * sin(x[0]) ) + 0.3
dL  = (cos(beta_1))/(sin(beta_1)**2)*db1*(l4 -l1*sin(x[0]) - 2 * r * cos(beta_1)) + (cos(beta_2))/(sin(beta_2)**2)*db2*(l4 -l1*sin(x[0]) - 2 * r * cos(beta_2)) - l1 * x[1] * cos(x[0]) * (csc(beta_1) + csc(beta_2))
L_o  = L / ( x[2] / (Y * A) + 1)
dL_o = r * (u[0] - u[1])


f[0] = x[1]
f[1] = (1/(I + l1 * m1 * r)) * (x[2] * l1 * (sin(x[0] + beta_1) - sin(x[0] - beta_2)) - l1 * 9.8 * (m1 + m2) * cos(x[0]) - k_f * x[1])
f[2] = (1/L_o)*(Y * A * dL - x[2] * dL_o - Y * A * dL_o) 

h = Matrix([x[0]])


#f1 = lambda x: scipy.tan(x) * (l1 * scipy.cos(float(X[0])) - l2 + l3 - 2 * r * scipy.sin(x) ) - (l4 - l1 * scipy.sin(float(X[0])) + 2 * r * scipy.cos(x))
#f2 = lambda x: scipy.tan(x) * (l2 - l1 * scipy.cos(float(X[0])) - 2 * r * scipy.sin(x) ) - (l4 - l1 * scipy.sin(float(X[0])) + 2 * r * scipy.cos(x))

#rb1 = scipy.optimize.newton(f1,scipy.optimize.brute(lambda x: f1(x)**2,((0,3.14),),Ns=100),tol=1e-20)[0]
#rb2 = scipy.optimize.newton(f2,scipy.optimize.brute(lambda x: f2(x)**2,((0,3.14),),Ns=100),tol=1e-20)[0]


#Kalman_Extendido(f,h,X,U,x,u,{diff(beta_1,x[0]):db1,diff(beta_2,x[0]):db2},{Function("beta1")(X[0]):rb1,Function("beta2")(X[0]):rb2,u[0]:U[0],u[1]:U[1]})
#Alta_Ganancia(f,h,X,U,x,u,{diff(beta_1,x[0]):db1,diff(beta_2,x[0]):db2},{})
Alta_Ganancia(f,h,X,U,x,u,{},{})


I   = 0.023
pprint(calcular_H([-2+4*I,-2-4*I,-10]))
pprint(calcular_H([-20+40*I,-20-40*I,-100]))
pprint(calcular_H([-200+400*I,-200-400*I,-1000]))
