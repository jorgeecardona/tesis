from sympy import *
from linear import Symbolic_Matrix

class hg:
    def parameters(self, s):
        return map(lambda x: var("(%s)"%(x)), s.split(" "))

    def state(self, n):
        self.x = Symbolic_Matrix("X_e", (n,1))
        return self.x

    def input(self, n):
        self.u = Symbolic_Matrix("U", (n,1))
        return self.u

    def f(self, f):
        self.fun_f = f

    def g(self, g):
        self.fun_g = g

    def h(self, h):
        self.fun_h = h

    functions = []

    def function(self, F, f, g, args):
        F = Function(F)
        self.functions.append((F, f, g, args))
        return F

    def hg(self):
        
        lie = self.fun_h

        M = zeros((self.fun_f.lines,1))
        for i in xrange(self.fun_f.lines):
            M[i] = lie[0]
            lie = lie.jacobian(self.x) * self.fun_f
            

        A = zeros(self.fun_f.lines)
        map(lambda x: A.mat.__setitem__(x,1), 
            range(1,self.fun_f.lines ** 2, self.fun_f.lines + 1))

        C = zeros((1,self.fun_f.lines))
        C[0] = 1
   
        # Jacobiano
        J = M.jacobian(self.x)

        for fun in self.functions:
            for args in fun[3]:
                for i in xrange(len(args)):

                    f = fun[1]
                    g = fun[2]
                    y = var("y")
                    F = (g(y,*args) - f(y,*args)).diff(args[i]) / (f(y,*args) - g(y,*args)).diff(y)

                    F = F.subs(y, fun[0](*args))
                    sub = {fun[0](*args).diff(args[i]): F}

                    J = J.subs(sub)

        J_inv = J.inv()
        
        return (M, J, J_inv)

    def compute_h(self, raices):
        n = len(raices)
        H = Symbolic_Matrix("H",(n,1))
        
    
        A = zeros(n)
        A[:-1,1:] = eye(n-1)
        C = zeros((1, n))
        C[0] = 1
        
        # Crear el polinomio a partir del determninante.
        var('lam')
        p = Poly((lam*eye(n) - A + H * C).berkowitz_det(),lam)

        # Crear el otro polinomio a partir de las raices.
        p1 = 1
        for i in raices:
            p1 *= Poly(lam-i, lam)

        # Comparar grados
        if p.degree != p1.degree:
            return None

        # Crear sistema de ecuaciones
        sys = ()
        for i in xrange(p.degree + 1):
            sys += ((p.coeff(i)-p1.coeff(i)),)

        r = solve(sys, list(H), simplified=False)
        
        # Reconstruir H con soluciones
        for i in xrange(H.lines):
            H[i,0] = r[H[i,0]]
            
        # Retornar H
        return H
        
