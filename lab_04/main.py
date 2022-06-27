import math
from numpy import e
import matplotlib.pyplot as plt
import numpy as np
from math import fabs
from mpl_toolkits.mplot3d import Axes3D

l = 10
T0 = 300
R = 0.5
k0 = 1
a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1

a2 = 2.049
b2 = 0.563e-3
c2 = 0.528e5
m2 = 1
alpha0 = 0.05
alphaN = 0.01

Fmax = 50
tmax = 60

EPS = 0.01
def lambd(T):
    return a1 * (b1 + c1 * pow(T, m1))

def c(T):
    return a2 + b2 * pow(T, m2) - c2 / (T * T)

def alpha(x):
    return 0.05 * l / (4 * x + l)

def F0(t):
    return (t / tmax) * Fmax * math.exp(-(t / tmax - 1))

def k(T):
    return lambd(T)

def p(x):
    return 2. / R * alpha(x)

def f(x, t, T):
    tmp = T / 300
    tmp2 = k0 * tmp * tmp
    return 2 * T0 / R * alpha(x) + tmp2 * F0(t) * math.exp(-tmp2 * x)

def solve_matrix(coefs, res):
    ksi = []
    eta = []

    n = len(coefs) - 1

    ksi.append(0)
    eta.append(0)

    ksi.append(-coefs[0][1] / coefs[0][0])
    eta.append(res[0] / coefs[0][0])


    for i  in range(1, n):
        a = coefs[i][0]
        b = -coefs[i][1]
        c = coefs[i][2]

        ksi.append(c / (b - a * ksi[i]))
        eta.append((-res[i] + a * eta[i]) / (b - a * ksi[i]))


    r = []
    r.append((res[n] - coefs[n][1] * eta[n]) / (coefs[n][1] * ksi[n] + coefs[n][2]))
    for i in range(n):
        r.append(r[i] * ksi[n - i] + eta[n - i]);

    for i in range(len(r) // 2):
        temp = r[i]
        r[i] = r[len(r) - i - 1]
        r[len(r) - i - 1] = temp

    return r
def check(a, b):
    for i in range(len(a)):
        if abs(a[i] - b[i] > EPS):
            return False
    return True

def boundary_condition(x0, xn, xstep):
    res = []
    x = x0
    while x < xn + EPS:
        res.append(T0)
        x += xstep
    return res
def dif_scheme(prev, prev_layer, t,  x0, xn, xstep, t0, tn, tstep):
    tres = []
    cres = []
    h = xstep
    tau = tstep

    tres.append( (
        (c(prev[0]) + (c(prev[0]) + c(prev[1])) / 4.) * h / 4. - alpha0 * tau + ((p(0) + p(h)) / 4. + p(0)) * tau * h / 4. + (k(prev[0]) + k(prev[1])) / (2. * h) * tau,
            (c(prev[0]) + c(prev[1])) * h / 16. - (k(prev[0]) + k(prev[1])) / (2. * h) * tau + (p(prev[0]) + p(prev[1])) * tau * h / 16.,
            0.))

    cres.append(((c(prev[0]) * prev_layer[0] + (c(prev[0]) + c(prev[1])) * (prev_layer[0] + prev_layer[1]) / 4.) * h / 4.
                    + alpha0 * T0 * tau +
                    (3. * f(0, t, prev[0]) + f(h, t, prev[1])) * h * tau / 8.))

    i = 1
    n = len(prev) - 2
    x = x0 + xstep
    for i in range(n):
        A = (k(prev[i - 1]) + k(prev[i])) * tau / h / 2.
        D = (k(prev[i + 1]) + k(prev[i])) * tau / h / 2.
        tres.append((
                A,
                -(A + D + c(prev[i]) * h + p(x) * tau * h),
                D
        ))

        cres.append(-(f(x, t, prev[i]) * tau * h + c(prev[i]) * prev_layer[i] * h))
        i += 1
        x += xstep


    n = len(prev) - 1

    tres.append((
            0.,
            ((c(prev[n - 1]) + c(prev[n])) / 4. + c(prev[n])) * h / 4. + (k(prev[n]) + k(prev[n - 1])) * tau / h / 2 + alphaN * tau + p(l) * tau * h / 4. + (p(l) + p(l - h)) * tau * h / 16.,
            (c(prev[n-1]) + c(prev[n])) * h / 16. - (k(prev[n]) + k(prev[n-1])) * tau / h / 2 + (p(l) + p(l-h)) * h * tau / 16.
        ))

    cres.append((((c(prev[n - 1]) + c(prev[n])) * (prev_layer[n-1] + prev_layer[n]) / 4. + c(prev[n]) * prev_layer[n]) * h / 4. + alphaN * T0 * tau + (3. * f(l, t, prev[n]) + f(l - h, t, prev[n-1])) * h * tau / 8.))
    return (tres, cres)

def iteration(prev_layer, t, x0, xn, xstep, t0, tn, tstep):
    prev = prev_layer
    coefs = dif_scheme(prev_layer, prev_layer, t, x0, xn, xstep, t0, tn, tstep)
    cur = solve_matrix(coefs[0], coefs[1])
    while !check(cur, prev):
        coefs = dif_scheme(prev, prev_layer, t, x0, xn, xstep, t0, tn, tstep)
        prev = cur
        cur = solve_matrix(coefs[0], coefs[1])
    return cur

def create_mesh(x0, xn, xstep, t0, tn, tstep):
    res = []
    res.append(boundary_condition(x0, xn, xstep))
    i = 0
    t = t0 + tstep
    while t < tn + EPS:
        res.append(iteration(res[i], t,  x0, xn, xstep, t0, tn, tstep))
        i += 1
        t += tstep
    return res


def main():
    x0 = 0
    xn = 10
    xstep = 1e-2
    t0 = 0
    tn = 10
    tstep = 1e-3
    
    res = create_mesh(x0, xn, xstep, t0, tn, tstep)
    
    f = open("out.txt", "w")
    f.writelines('\n'.join(str(x) for x in res))
    f.close()
    
    x = []
    xi = x0
    while xi < xn + EPS:
        x.append(xi)
        xi += xstep

    te = []
    
    t = t0
    while t < tn + EPS:
        te.append(t)
        t += tstep

    arrs = np.array(res)


    xgrid, ygrid = np.meshgrid(x, te)
    xgrid = np.array(xgrid)
    ygrid = np.array(ygrid)
    zgrid = np.array(arrs)

    fig = plt.figure(figsize=(7, 4))
    ax_3d = Axes3D(fig)
    ax_3d.plot_surface(ygrid, xgrid, zgrid, cmap='GnBu')
    ax_3d.set_ylabel("x, cm")
    ax_3d.set_xlabel("t, c")
    ax_3d.set_zlabel("T, K")
    plt.show()
main()
