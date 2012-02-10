# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or    
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of    
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License    
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from sympy import Symbol, var, Function, zeros, nsolve

def Symbolic_Matrix(name,size):
    M = zeros(size)
    for i in xrange(M.lines):
        for j in xrange(M.cols):
            M[i,j] = Symbol("%s[%i,%i]"%(name,i,j))
    return M

class linear:

    def parameters(self, s):
        return map(lambda x: var("p.%s"%(x)), s.split(" "))

    def state(self, n):
        self.x = Symbolic_Matrix("X_e", (n,1))
        return self.x

    def input(self, n):
        self.u = Symbolic_Matrix("U", (n,1))
        return self.u

    def f(self, f):
        self.fun_f = f
    
    def h(self, h):
        self.fun_h = h

    functions = []

    def function(self, F, f, g, args):
        F = Function(F)
        self.functions.append((F, f, g, args))
        return F

    def linearize(self, x_o, u_o):
        if x_o.lines != self.x.lines:
            raise Exception, "Different line size"

        if x_o.cols != self.x.cols:
            raise Exception, "Different line size"

        if u_o.lines != self.u.lines:
            raise Exception, "Different line size"

        if u_o.cols != self.u.cols:
            raise Exception, "Different line size"


        A = self.fun_f.jacobian(self.x)
        B = self.fun_f.jacobian(self.u)
        C = self.fun_h.jacobian(self.x)


        for F in self.functions:
            

            sub ={F[0]}
        
        for fun in self.functions:
            for args in fun[3]:
                for i in xrange(len(args)):

                    f = fun[1]
                    g = fun[2]
                    y = var("y")
                    x = [var("x_%d"%(k) for k in xrange())]



                    F = (
                        g(y,*args) - f(y,*args)).diff(args[i]) / (
                        f(y,*args) - g(y,*args)).diff(y)

                    F = F.subs(y, fun[0](*args))
                    sub = {fun[0](*args).diff(args[i]): F}

                    A = A.subs(sub)
                    B = B.subs(sub)
                    C = C.subs(sub)


        num = True


        for fun in self.functions:
            for args in fun[3]:

                sub = dict(map(
                        lambda i: (self.x[i],x_o[i]),
                        xrange(self.x.lines * self.x.cols))
                           + map(
                        lambda i: (self.u[i],u_o[i]),
                        xrange(self.u.lines * self.u.cols)))

                args_o = map(lambda x: x.subs(sub),args)

                num = True
                for arg in args_o:
                    num &= arg.is_number

                if num:
                                        
                    y = var("y")
                    y = nsolve(
                        fun[1](y,*args_o)- fun[2](y,*args_o),
                        0.7)

                    A = A.subs(fun[0](*args), y)
                    B = B.subs(fun[0](*args), y)
                    C = C.subs(fun[0](*args), y)


        sub = dict(map(
                lambda i: (self.x[i],x_o[i]),
                xrange(self.x.lines * self.x.cols))
                   + map(
                lambda i: (self.u[i],u_o[i]),
                xrange(self.u.lines * self.u.cols)))
        
        A = A.subs(sub).evalf()
        B = B.subs(sub).evalf()
        C = C.subs(sub).evalf()

        return  (A,B,C)
