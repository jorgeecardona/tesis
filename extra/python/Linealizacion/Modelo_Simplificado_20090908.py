

from Linealizacion import *
from sympy import Matrix, Function, cos, sin, tan, diff, pprint
import scipy
import scipy.optimize
sec = lambda x: 1 / cos(x)


x = Symbolic_Matrix("x",(3,1))
u = Symbolic_Matrix("u",(2,1))

X = Matrix([0 , 0 , 0.05166])
U = Matrix([2, 2.307515])

M_p     = 0.86
Y       = 2e6
A       = 14e-6
l       = 0.41
r       = 0.05
k_f     = 1
L_k     = 0.3
rho_o   = 714
rho_e   = 714

beta, L, T, rho = map(Function, ["beta", "L", "T", "rho"])

L = lambda beta: L_k + 4 * r * ( beta - tan(beta)) + l * sec(beta)
T = lambda rho: Y * A * (rho_o /rho - 1)
rho = lambda L, M: M /(A * L)

f = zeros((3,1))

f[0] = x[1]
f[1] = 9.8 - 2 * T(rho(L(beta(x[0])),x[2])) * sin(beta(x[0]))/M_p - k_f * x[1] / M_p
f[2] = A * r * (rho_e * u[0] - rho(L(beta(x[0])),x[2]) * u[1] )

h = Matrix([x[0]])

db1 = 1 / ((l / 2) * sec(beta(x[0]))**2 - 2 * r * tan (beta(x[0])) * sec(beta(x[0])))
Lineal(f,h,X,U,x,u,{diff(beta(x[0]),x[0]):db1},{})



f1 = lambda x: l * scipy.sin(x) - 4 * r  - 2 * float(X[0]) * scipy.cos(x)
rb1 = scipy.optimize.newton(f1,scipy.optimize.brute(lambda x: f1(x)**2,((0,3.14),),Ns=100),tol=1e-20)[0]


db1 = 1 / ((l / 2) * sec(beta(x[0]))**2 - 2 * r * tan (beta(x[0])) * sec(beta(x[0])))
Lineal(f,h,X,U,x,u,{diff(beta(x[0]),x[0]):db1},{beta(X[0]):rb1})

#print f.jacobian(x).subs(beta,beta(x[0]))
#print f.jacobian(u)
#print h.jacobian(x)


