import numpy as np
from scipy.integrate import quadrature

def A(N):
    return np.vander(np.linspace(-1,1,N), increasing=True)

def B(N):
    return np.eye(N)

def coefs(N):
    return np.linalg.solve(A(N), B(N))

def poly(ii, N):
    alpha = coefs(N)
    return np.poly1d(alpha.T[ii][::-1])

def polys(N):
    p = []
    for i in range(N):
        p.append(poly(i, N))
    return p

def polyder(ii, N):
    alpha = coefs(N)
    p = np.poly1d(alpha.T[ii][::-1])
    return np.polyder(p)

def polyders(N):
    p = []
    for i in range(N):
        p.append(polyder(i, N))
    return p

def getX(coords):
    N = len(coords)
    def x(s):
        y = 0
        p = polys(N)
        for ii in range(N):
            y += p[ii](s)*coords[ii]
        return y
    return x

def getJ(coords):
    N = len(coords)
    def j(s):
        y = 0
        p = polyders(N)
        for ii in range(N):
            y += p[ii](s)*coords[ii]
        return y
    return j

def getIe(coords):
    N = len(coords)
    X = getX(coords)
    J = getJ(coords)
    p0 = polys(N)
    p1 = polyders(N)

    Ie = np.empty((N, N))
    for ii in range(N):
        for jj in range(N):
            fun1 = lambda s : p1[ii](s)/J(s)
            fun2 = lambda s : p1[jj](s)/J(s)
            fun  = lambda s : fun1(s)*fun2(s)*J(s)
            Ie[ii, jj] = quadrature(fun, -1, 1)[0]

    return Ie
