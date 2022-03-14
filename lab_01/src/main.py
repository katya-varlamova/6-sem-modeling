import matplotlib.pyplot as plt
#pikar
def int_1(x, yn = 0, h = 0):
    x = x + h
    return x * x * x / 3
def int_2(x, yn = 0, h = 0):
    x = x + h
    return int_1(x) + pow(x, 7) / 63
def int_3(x, yn = 0, h = 0):
    x = x + h
    return int_2(x) + 2 * pow(x, 11) / 2079 + pow(x, 15) / 59535
def int_4(x, yn = 0, h = 0):
    x = x + h
    return int_3(x) + 13 * pow(x, 15) / 218295 + 82 * pow(x, 19) / 37328445 + 662 * pow(x, 23) / 10438212015 + 4 * pow(x, 27) / 3341878155 + pow(x, 31) / 109876902975;
# runge
def fi(x, y):
    return x * x + y * y
alpha = 0
def runge(xn, yn, h):
    f = fi(xn, yn)
    tmp =  h / (2 * alpha)
    return yn + h * ((1 - alpha) * f + alpha * fi(xn + tmp, yn + tmp * f))

# euler
def eul(xn, yn, h):
    return yn + h * fi(xn, yn)
st = 5
def calc_results(x0, b, h, y0, funcs):
    results = []
    args = []
    x = x0
    mul = pow(10, st)
    mod = pow(10, st - 2)
    while x < b:
        if int((x * mul) % mod) == 0:
            args.append(x)
        x += h
    results.append(args)
    for f in funcs:
        tmp = []
        yl = y0
        x = x0
        while x < b:
            if int((x * mul) % mod) == 0:
                tmp.append(yl)
            yl = f(x, yl, h)
            x += h
        results.append(tmp)
    return results

def print_table(r):
    print("x" + "-" * 55 + "x")
    print("|{:^7}|{:^7}|{:^7}|{:^7}|{:^7}|{:^7}|{:^7}|".format("x", "p1", "p2", "p3", "p4", "run", "eul"))
    print("x" + "-" * 55 + "x")
    for i in range(len(r[0])):
        print("|{:^7.2f}|{:^7.2f}|{:^7.2f}|{:^7.2f}|{:^7.2f}|{:^7.2f}|{:^7.2f}|".format(r[0][i], r[1][i], r[2][i], r[3][i], r[4][i], r[5][i], r[6][i]))
    print("x" + "-" * 55 + "x")
def view_graph(r):
    for i in range(len(r)):
        neg = [-i for i in r[i][:0:-1]]
        r[i] = neg + r[i]
    plt.xlabel("X")
    plt.ylabel("Y")
    names = ["pikar1", "pikar2", "pikar3", "pikar4", "runge", "euler"]
    for i in range(1, len(r)):
        plt.plot(r[0], r[i], label = names[i - 1])
    plt.legend()
    plt.show()
alpha = 0.5
r = calc_results(0.0, 2.2, pow(10, -st), 0, [int_1, int_2, int_3, int_4, runge, eul])
print_table(r)
view_graph(r)
