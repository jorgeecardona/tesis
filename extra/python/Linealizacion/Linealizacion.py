
from sympy import var, diff, Matrix, zeros, sin, pprint,tan, latex


def Symbolic_Matrix(name,size):
    M = zeros(size)
    for i in xrange(M.lines):
        for j in xrange(M.cols):
            M[i,j] = var("%s_%i_%i"%(name,i,j))
    return M



def Substitute_Matrix(*args,**kwords):
    expression = args[0]
    ch = {}
    for k in xrange(1,len(args),2):
        for i in xrange(args[k].lines):
            for j in xrange(args[k].cols):
                ch[args[k][i,j]] = args[k+1][i,j]
    return expression.subs(ch)




def Lineal(f,h,x,u,var_x,var_u,d,v):
    A = zeros((3,3))
    A[0]

    # Calcular matrices
    A = Substitute_Matrix(f.jacobian(var_x).subs(d),var_x,x,var_u,u).subs(v)
    B = Substitute_Matrix(f.jacobian(var_u).subs(d),var_x,x,var_u,u).subs(v)
    C = Substitute_Matrix(h.jacobian(var_x).subs(d),var_x,x,var_u,u).subs(v)


    pprint(A)
    pprint(B)
    pprint(C)

