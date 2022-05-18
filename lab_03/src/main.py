import matplotlib.pyplot as plt
import math
import pylab
from matplotlib import ticker
R0 = 0.35
Tw = 2000
T0 = 10000
k0 = 0.1
c = 3e10
p = 4
h = 0.001

def T(z):
    return (Tw - T0) * pow(z, p) + T0

def k(z):
    tmp = T(z) / 300
    return k0 * tmp * tmp

def up(z):
    return (3.084 * 1e-4) / (math.exp(47990 / T(z)) - 1)

def K(z):
    return c / (3 * k(z))

def P(z):
    return c * k(z)

def F(z):
    return -c * k(z) * up(z)

def X(z, h):
    return (K(z) + K(z + h)) / 2

def V(z, h):
    tmp1 = (z - h / 2) * (z - h / 2)
    tmp2 = (z + h / 2) * (z + h / 2)
    return (tmp2 - tmp1) / 2


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

def dif_scheme(z, h, stop):
    coefs = []
    res = []
    coefs.append([-1. / R0 / R0 / h * X(z, h) * (z + h/2.) - h / 4. * V(z + h / 4., h / 2.) * (P(z) + 1 / 2. * P(z + h / 2)),
                  1. / R0 / R0 / h * X(z, h) * (z + h/2.) - h / 8. * V(z + h/4., h/2.) * P(z + h / 2.),
                  0])
    res.append(((F(z) + F(z + h/2.)) * h * V(z + h / 4., h / 2.) / 4.))

    z += h
    
    while z < stop - h / 2:
        A = (z - h / 2.) * X(z - h, h) / (R0 * R0 * h)
        C = (z + h / 2.) * X(z, h) / (R0 * R0 * h)

        coefs.append([A, -(A + C + P(z) * V(z, h)), C])
        res.append(F(z) * V(z, h))
        z += h

    coefs.append([0,
                  (z - h / 2.) * X(z - h, h) / (R0 * R0 * h) - h * V(z - h / 4., h / 2.) * P(z - h / 2.) / 8.,
                  -(z - h / 2.) * X(z - h, h) / (R0 * R0 * h) - z * 0.393 * c / R0 - h / 4. * V(z - h / 4., h / 2.) * (P(z - h/2.) / 2. + P(z))
                  ])
    res.append(h * (F(z - h / 2.) + F(z)) * V(z - h / 4., h / 2.))

    return (coefs, res)


def c_F_simpson(z, zi, u):
    return k(z) * (up(z) - u[zi]) * z
def simpson(func, a, b, num_of_nodes, u):
    if num_of_nodes < 3:
        num_of_nodes = 3
    if num_of_nodes & 1 == 0:
        num_of_nodes -= 1

    h = (b - a) / (num_of_nodes - 1)
    x = a
    res = 0

    for i in range((num_of_nodes - 1) // 2):
        res += func(x, i, u) + 4 * func(x + h, i + 1, u) + func(x + 2 * h, i + 2, u)
        x += 2 * h

    return res * (h / 3)
def def_F_simpson(z, u):
    f = []
    for i in range(len(z)):
        try:
            f.append(c * R0 / z[i] * simpson(c_F_simpson, 0, z[i], i, u))
        except ZeroDivisionError:
            f.append(0)
    return f
                 
def center(y, h, i):
    return (y[i + 1] - y[i - 1]) / 2 / h
def c_f(x, der):
    return -K(x) * der / R0
def def_F(x, u, h):
    f = [0]
    n = len(u)
    for i in range(1, n - 1):
        f.append(c_f(x[i], center(u, h, i)))
    f.append(c_f(x[n - 1], (-4 * u[n - 2] + u[n - 3] + 3 * u[n - 1]) / 2 / h))
    return f


def solve(arg):
    r = dif_scheme(0, h, arg)
    u = solve_matrix(r[0], r[1])
    return u[len(u) - 1]

def second_der(y, h, i):
    return (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h) if i < len(y) - 1 and i > 0 else "-"
def divF(x, u, h):
    f = [0]
    n = len(u)
    for i in range(1, n - 1):
        f.append(c_f(x[i], second_der(u, h, i)))
    return f

def trace_R0(init, end):
    step = (end - init) / 100
    global R0
    R0 = init
    val = []
    us = []
    for _ in range (100):
        val.append(R0)
        us.append(solve(1))
        R0 += step
    return (val, us)

def trace_Tw(init, end):
    step = 500
    global Tw
    Tw = init
    val = []
    us = []
    for _ in range (8):
        val.append(Tw)
        us.append(solve(1))
        Tw += step
    return (val, us)

def trace_T0(init, end):
    step = (end - init) / 10
    global T0
    T0 = init
    val = []
    us = []
    for _ in range (10):
        val.append(T0)
        us.append(solve(1))
        T0 += step
    return (val, us)

def trace_k0(init, end):
    step = (end - init) / 10
    global k0
    k0 = init
    val = []
    us = []
    for _ in range (10):
        val.append(k0)
        us.append(solve(1))
        k0 += step
    return (val, us)

def trace_p(init, end):
    step = (end - init) / 11
    global p
    p = init
    val = []
    us = []
    for _ in range (11):
        val.append(p)
        us.append(solve(1))
        p += step
    return (val, us)
def main():
    r = dif_scheme(0, h, 1)
    u = solve_matrix(r[0], r[1])

    z = []
    z_i = 0
    up_v = []
    for i in range(len(u)):
        z.append(z_i)
        up_v.append(up(z_i))
        z_i += h
        
    f1 = def_F(z, u, h)
    f2 = def_F_simpson(z, u)

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1, 1)) 

    pylab.figure(1)
    ax1 = pylab.subplot(1, 3, 1)
    pylab.plot(z, f1, label = "f(z) - numerical differentiation")
    pylab.legend()
    pylab.plot(z, f2, label = "f(z) - integration (simpson method)")
    pylab.legend()

    ax2 = pylab.subplot(1, 3, 2)
    pylab.plot(z, u, label = "u(z)")
    pylab.plot(z, up_v, label = "up(z)")
    pylab.legend()

    ddf = divF(z, u, h)
    ax3 = pylab.subplot(1, 3, 3)
    pylab.plot(z[1:-1], ddf[1:], label = "divF(z)")
    pylab.legend()

    ax2.yaxis.set_major_formatter(formatter)
    plt.sca(ax1)
    ax1.yaxis.set_major_formatter(formatter)
    plt.sca(ax3)
    ax3.yaxis.set_major_formatter(formatter)
    pylab.show()

    pylab.figure(2)

    val, us = trace_R0(1, 101)
    global R0
    R0 = 0.35
    pylab.subplot(2, 3, 1)
    pylab.plot(us, val, label = "u0(R)")
    pylab.legend()

    val, us = trace_Tw(500, 4000)
    global Tw
    Tw = 2e3
    pylab.subplot(2, 3, 2)
    pylab.plot(val,us, label = "u0(Tw)")
    pylab.legend()

    val, us = trace_k0(0.0004, 0.0014)
    global k0
    k0 = 8e-4
    pylab.subplot(2, 3, 3)
    pylab.plot(val, us, label = "u0(k0)")
    pylab.legend()

    val, us = trace_T0(5000, 15000)
    global T0
    T0 = 1e4
    pylab.subplot(2, 3, 4)
    pylab.plot(val, us, label = "u0(T0)")
    pylab.legend()

    val, us = trace_p(4, 15)
    global p
    p = 4
    pylab.subplot(2, 3, 5)
    pylab.plot(val, us, label = "u0(p)")
    pylab.legend()

    pylab.show()

main()
