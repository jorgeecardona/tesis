from numpy import array, matrix, sin, sqrt, dot, cos, ix_, zeros, concatenate, abs, log10, exp, ones
from numpy.linalg import norm

from mpmath import mpf, mp
mp.dps=200

from normal import normal

def minimum(P):
#    print "Calling find_min with P: ", P

    if len(P) == 1:
        return P[0]

    eps = mpf(10)**-12

    P = [array([mpf(i) for i in p]) for p in P]
    
    # Step 0. Choose a point from C(P)
    x  = P[array([dot(p,p) for p in P]).argmin()]

    for __i__ in xrange(1000):
#        print x

        # Step 1. \alpha_k := min{x_{k-1}^T p | p \in P}
        p_alpha = P[array([dot(x,p) for p in P]).argmin()]

        if dot(x,x-p_alpha) < eps:
            return array([float(i) for i in x]) 
        
        Pk = [p for p in P if abs(dot(x,p-p_alpha)) < eps]

        # Step 2. P_k := { p | p \in P and x_{k-1}^T p = \alpha_k}
        P_Pk = [p for p in P if not array([(p == q).all() for q in Pk]).any()]

        if len(Pk) == len(P):
            return array([float(i) for i in x]) 

        y = minimum(Pk)


        p_beta = P_Pk[array([dot(y,p) for p in P_Pk]).argmin()]
        
        if dot(y,y-p_beta) < eps:
            return array([float(i) for i in y]) 

        # Step 4. 
        P_aux = [p for p in P_Pk if (dot(y-x,y-p)>eps) and (dot(x,y-p)!=0)]
        p_lambda = P_aux[array([dot(y,y-p)/dot(x,y-p) for p in P_aux]).argmin()]
        lam = dot(x,p_lambda-y) / dot(y-x,y-p_lambda)

        if lam < eps:
            return array([float(i) for i in x]) 

        x += lam * (y-x)

    print "I'm sick of this error."
    return (x + y) / 2

dynmap={}
def dynamic_programming(func):
    def decorator(P):
        k = tuple(tuple(p) for p in P)
        if k in dynmap:
            return dynmap[k]
        else:
            r = func(P)
            dynmap[k] = r
            return r 
    return decorator


@dynamic_programming
def __minimum(P, eps = 1e-8):

    # Number of elements.
    n = len(P)
    if n == 0:
        raise Exception, "0 points"
    
    if n == 1:
        return  P[0]

    m = len(P[0])

    # Same size
    if array([len(p) != m for p in P]).any():
        raise Exception, "Inconsistent point sizes."

    # A possible face
    if n == m:
        return wolfe(P)

    if n < m:
        return wolfe(P)
    
    mins = []

    for i in xrange(len(P)):
#        print "P_%d: "%(i), P[i]

        def combine(P, n):

            if len(P) < n:
                return []

            if n == 1:
                return [[p] for p in P]
            
            if n == len(P):
                return [P]

            return [[P[0]] + p for p in combine(P[1:], n-1)] + combine(P[1:], n)

        # Select n-1 points to create a face candidate.
        Pk_1 = P[:i] + P[i+1:]
        Pk = [[P[i]] + p for p in combine(Pk_1, m-1)]
        
        # Pk is a list of face candidates.
#        print P
#        print Pk

        faces = []
        Pk_2 = [P[i]]
        for p in Pk:
            dot_P = array([dot(normal(p),q - p[0]) for q in Pk_1])
            # print "p: ", p, normal(p), dot_P

            if dot_P.min() * dot_P.max() >= -eps:
                faces.append(p)

                [Pk_2.append(q) for q in p[1:] if not array([(q == r).all() for r in Pk_2]).any()]
                # print "Cara: ", p
                # print dot_P.min()
                # print dot_P.max()
        
        # print "len(P): ", len(P)
        # print "len(Pk_2): ", len(Pk_2)

        mins.append(wolfe(Pk_2))
        mins.append(wolfe(Pk_1))


    return  mins[array([dot(p,p) for p in mins]).argmin()]    

wolfe = minimum
minimum = __minimum

if __name__ == '__main__':

    print __minimum([array([  5.51404676,  -9.66783678,   7.22818865,  18.25147299,
        19.94377563, -16.29059911]), array([ -7.59458805,   5.15059326, -11.15044168,   3.12526983,
         6.85910303,   1.97949357]), array([  8.7489085 , -31.53282378,  -4.22016299,  48.99492468,
       -25.63744762, -11.63770461]), array([ 37.50040804,  15.64426383, -11.81390286, -16.81087638,
        29.06938081,  -4.30141625]), array([ 20.35506981,   6.5493941 ,   5.71456202,  29.0666551 ,
       -21.9927714 ,   7.09042777]), array([ -2.71856809,  10.95897922,  20.11699329,   9.87328894,
         4.01778286, -39.65485385]), array([ 10.25613042,  -4.87234234,  12.20360493,  18.53197223,
         7.71111023,  28.53474336]), array([ 11.91368994,  33.33886916, -73.75246764, -39.8114611 ,
        -9.29275498,  -9.68341651]), array([  8.39065272, -12.64968497, -20.15945961, -21.80692032,
        25.84139455, -23.56411231]), array([ 19.71619665,  15.90476122, -29.97520753, -66.38031366,
        29.01398409, -36.36284803]), array([ 14.10557533, -22.14770093, -11.48532974,  22.08802499,
       -30.58627905, -24.72641955])])

    print __minimum([
            array([-3.0668369 ,  0.99226292,  1.42386285,  0.92769034, -0.13652053, 0.34232693]), 
            array([ 1.5830313 ,  0.36909878, -4.65754129, -1.0066756 ,  0.44759627, 0.39416678]), 
            array([-0.26126152, -0.81416674, -0.09124101, -0.06502509, -1.23242146, -0.6132118 ]), 
            array([-0.26225631,  0.97401745, -0.94403595,  0.25495119,  1.8674588 , 1.85616728]), 
            array([ 0.24847795,  0.02949977, -0.53689046,  0.06477375, -0.96588737, -0.91785888]), 
            array([ -1.28394310e+00,  -2.10186025e-01,   1.06033684e+00, -1.14145396e-01,   1.24223681e-03,   5.17836481e-01]), 
            array([ 0.99138383, -0.06205938, -1.40463356, -0.31555266,  0.66419493, 0.72847102])
            ])

    print __minimum([array([ 1.05880179, -0.31739554, -1.3659017 ,  0.28531694,  0.80913271,
        0.44038549]), array([-3.80494612, -0.25725629, -0.65919596,  0.16979377, -2.06476782,
       -1.25979309]), array([-1.80453767, -1.14825985,  1.0889194 , -0.01448414, -0.04681101,
       -0.75456722]), array([ 0.72129753,  0.01714778, -0.49100158, -0.04538874, -0.16723215,
        0.27349124]), array([-0.77345019, -0.38400753,  0.863607  ,  0.2595667 ,  1.57016811,
        1.35950704]), array([-0.80466675,  0.22005762,  0.4010332 ,  0.25371974,  0.03755512,
        0.04845693]), array([ 0.87044265, -0.20873802, -1.87822948, -0.03636987, -1.67328671,
       -1.6663354 ]), array([ 1.44856219,  0.17636973, -1.31034667, -0.66606963,  0.0362734 ,
        0.25616458])])


    print __minimum([array([1,1]), array([3,0]), array([0,3]), array([3,3]), array([4,4])])
    print __minimum( [array([ -4.83907292e+00,   2.22438863e+04,  -2.67496763e+04]), array([   9.71147604, -351.46404195, -292.18064276]), array([  4.60452808e+00,   1.07020174e+05,  -1.25310230e+05]), array([  2.16080134e+00,   5.12019937e+04,  -5.96167833e+04]), array([  2.65472146e+00,   6.70546443e+04,  -7.71619656e+04]), array([  1.55775358e+00,  -1.34347516e+05,   1.53209265e+05]), array([   13.22464295,  1869.01251292, -2137.61850989])])


    print __minimum( [array([ -4.83907292e+00,   2.22438863e+04,  -2.67496763e+04]), array([   9.71147604, -351.46404195, -292.18064276]), array([  4.60452808e+00,   1.07020174e+05,  -1.25310230e+05]), array([  2.16080134e+00,   5.12019937e+04,  -5.96167833e+04]), array([  2.65472146e+00,   6.70546443e+04,  -7.71619656e+04]), array([  1.55775358e+00,  -1.34347516e+05,   1.53209265e+05]), array([   13.22464295,  1869.01251292, -2137.61850989]), array([ 12273.18670123,  -1233.32015854,  61690.10864825])])

    print __minimum([array([ -72024.29139685,     419.29373955,  110723.5960603 ]), array([    7.3549799 , -1255.92108206,     9.45781838]), array([ -1.85553318e+00,  -1.77497921e+02,  -1.14172014e+04])])

    print __minimum([array([  97316.2958198 ,     997.7025733 ,  151350.30532618]), array([ -6.25809145e-01,  -1.39333078e+03,   9.33530701e+00]), array([  6.26509323e+00,  -1.08298220e+03,  -8.11827490e+04]), array([ -1.29414885e+01,  -9.69105418e+02,   2.53632708e+05]), array([  7.82322896e+00,   1.88815305e+01,   5.48789950e+04]), array([ -1.28940611e+01,   5.17422301e+02,   6.38864539e+04]), array([ -8.16746785e+00,  -1.50394858e+03,  -1.36345350e+05]), array([  1.75609676e+01,  -8.91326957e+02,  -8.82899352e+04]), array([  7.39719564e+00,   1.74080382e+03,   1.88046341e+05]), array([  1.14821503e+01,   1.57353035e+03,   1.64662987e+05]), array([  2.19713781e+01,   3.02818173e+03,   3.09260435e+05]), array([  7.70859577e+00,   2.81759776e+03,   2.85362716e+05]), array([ -2824.69805732,  -2264.29013682,  61710.26692425]), array([  -5901.52425596,   -1729.72615394,  424085.64618553]), array([  -2699.36351202,    -712.81248933,  200460.47929018]), array([  3.40105331e+00,   6.45628057e+03,   6.48383854e+05])])

    print __minimum( [array([  97316.2958198 ,     997.7025733 ,  151350.30532618]), array([ -6.25809145e-01,  -1.39333078e+03,   9.33530701e+00]), array([  6.26509323e+00,  -1.08298220e+03,  -8.11827490e+04]), array([ -1.29414885e+01,  -9.69105418e+02,   2.53632708e+05]), array([  7.82322896e+00,   1.88815305e+01,   5.48789950e+04]), array([ -1.28940611e+01,   5.17422301e+02,   6.38864539e+04]), array([ -8.16746785e+00,  -1.50394858e+03,  -1.36345350e+05]), array([  1.75609676e+01,  -8.91326957e+02,  -8.82899352e+04]), array([  7.39719564e+00,   1.74080382e+03,   1.88046341e+05]), array([  1.14821503e+01,   1.57353035e+03,   1.64662987e+05]), array([  2.19713781e+01,   3.02818173e+03,   3.09260435e+05]), array([  7.70859577e+00,   2.81759776e+03,   2.85362716e+05]), array([ -2824.69805732,  -2264.29013682,  61710.26692425]), array([  -5901.52425596,   -1729.72615394,  424085.64618553]), array([  -2699.36351202,    -712.81248933,  200460.47929018]), array([  3.40105331e+00,   6.45628057e+03,   6.48383854e+05]), array([ -1.40618264e+01,  -2.78958249e+02,  -2.65327152e+04])])


