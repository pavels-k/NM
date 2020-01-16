import numpy as np
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math as mth
import matplotlib.pyplot as plt



x0 = 0
xl = np.pi / 4

y0 = 0
yl = np.log(2)


def progonka(a, b, c, d, s):
    P = np.zeros(s)
    Q = np.zeros(s)

    P[0] = -c[0] / b[0]
    Q[0] = d[0] / b[0]

    k = s - 1

    for i in range(1, s):
        P[i] = -c[i] / (b[i] + a[i] * P[i - 1])
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1])
    P[k] = 0
    Q[k] = (d[k] - a[k] * Q[k - 1]) / (b[k] + a[k] * P[k - 1])

    x = np.zeros(s)
    x[k] = Q[k]

    for i in range(s - 2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]
    return x


def Analitic(x, y, t, a):
    return np.cos(2 * x) * np.cosh(y) * np.exp(-3 * t * a)


def Fy0(x, t, a):
    return np.cos(2 * x) * np.exp(-3 * a * t)


def FyN(x, t, a):
    return (3 / 4) * np.cos(2 * x) * np.exp(-3 * a * t)


def Fx0(y, t, a):
    return np.cosh(y) * np.exp(-3 * a * t)


def FxN(y, t, a):
    return 0


def Psi(x, y):
    return np.cos(2 * x) * np.cosh(y)


def Alternating_Directions(hx, hy, r,  t, par_a):  # Переменных направлений
    x = np.arange(x0, xl + hx, hx)
    y = np.arange(y0, yl + hy, hy)
    time = np.arange(0, t + r, r)

    U = np.zeros((len(time), len(x), len(y)))

    for i in range(len(x)):
        for j in range(len(y)):
            U[0][i][j] = Psi(x[i], y[j]) #  5-ое условие
    for k in range(0, len(time)):
        for j in range(len(y)):
            U[k][0][j] = Fx0(y[j], time[k], par_a) # 1-ое условие 
            U[k][-1][j] = FxN(y[j], time[k], par_a) # 2-ое условие

    for k in range(0, len(time)):
        for i in range(len(x)):
            U[k][i][0] = Fy0(x[i], time[k], par_a) # 3-условие

    for k in range(1, len(time)):
        k_half_table = np.zeros((len(x), len(y)))

        for i in range(1, len(x) - 1):
            a = np.zeros(len(y))
            b = np.zeros(len(y))
            c = np.zeros(len(y))
            d = np.zeros(len(y))

            tmp = (par_a * r) / (hx**2**2 * 2)
            # первый дробный шаг
            for j in range(1, len(y) - 1):
                a[j] = tmp
                b[j] = -2 * tmp - 1
                c[j] = tmp

                # явная по всем значениям у
                d[j] = (- par_a * r / (hy**2 * 2)) * \
                    (U[k - 1][i][j + 1] - 
                     2 * U[k - 1][i][j] +
                     U[k - 1][i][j - 1]) -\
                    U[k - 1][i][j]

            alpha = 0
            betta = 1
            gamma = 1
            delta = 0
            b[0] = betta - alpha / hy
            c[0] = alpha / hy
            d[0] = Fy0(x[i], time[k] - r / 2, par_a)

            a[-1] = - gamma / hy
            b[-1] = delta + gamma / hy
            d[-1] = FyN(x[i], time[k] - r / 2, par_a)

            ans = progonka(a, b, c, d, len(d))
            k_half_table[i] = ans
            for j in range(len(y)):
                k_half_table[0][j] = Fx0(y[j], time[k] - r / 2, par_a)
                k_half_table[-1][j] = FxN(y[j], time[k] - r / 2, par_a)
        # второй дробный шаг
        for j in range(1, len(y) - 1):
            a = np.zeros(len(x))
            b = np.zeros(len(x))
            c = np.zeros(len(x))
            d = np.zeros(len(x))

            tmp = (par_a * r) / (hx**2 * 2)
            for i in range(1, len(x)):
                a[i] = tmp
                b[i] = -2 * tmp - 1
                c[i] = tmp
                # неявно по y
                d[i] = (- par_a * r / (hy**2 * 2)) * \
                    (k_half_table[i][j + 1] -
                     2 * k_half_table[i][j] +
                     k_half_table[i][j - 1]) -\
                    k_half_table[i][j]

            alpha = 0
            betta = 1
            gamma = 0
            delta = 1
            b[0] = betta - alpha / hx
            c[0] = alpha / hx
            d[0] = Fx0(y[j], time[k], par_a)

            a[-1] = - gamma / hx
            b[-1] = delta + gamma / hx
            d[-1] = FxN(y[j], time[k], par_a)

            ans = progonka(a, b, c, d, len(d))
            for i in range(len(ans)):
                U[k][i][j] = ans[i]
            for j in range(len(y)):
                U[k][0][j] = Fx0(y[j], time[k], par_a)
                U[k][-1][j] = FxN(y[j], time[k], par_a)

            for i in range(len(x)):
                U[k][i][0] = Fy0(x[i], time[k], par_a)
                U[k][i][-1] = FyN(x[i], time[k], par_a) * hy\
                    + U[k][i][-2] # производная

    for j in range(len(y)):
        U[-1][0][j] = Fx0(y[j], time[-1], par_a)
        U[-1][-1][j] = FxN(y[j], time[-1], par_a)

    for i in range(len(x)):
        U[-1][i][0] = Fy0(x[i], time[-1], par_a)

    return U[-1].transpose()


def Fractional_Steps(hx, hy, r,  t, par_a):
    x = np.arange(x0, xl + hx, hx)
    y = np.arange(y0, yl + hy, hy)
    time = np.arange(0, t + r / 2, r)

    U = np.zeros((len(time), len(x), len(y)))
    
    for i in range(len(x)):
        for j in range(len(y)):
            U[0][i][j] = Psi(x[i], y[j]) # 5

    for k in range(0, len(time)):
        for j in range(len(y)):
            U[k][0][j] = Fx0(y[j], time[k], par_a) # 1
            U[k][-1][j] = FxN(y[j], time[k], par_a) # 2

    for k in range(0, len(time)):
        for i in range(len(x)):
            U[k][i][0] = Fy0(x[i], time[k], par_a) # 3


    for k in range(1, len(time)):
        k_half_table = U[k].copy()
        for j in range(1, len(y) - 1):

            a = np.zeros(len(x))
            b = np.zeros(len(x))
            c = np.zeros(len(x))
            d = np.zeros(len(x))

            tmp = par_a * r / hx**2
            # первый дробный шаг
            for i in range(1, len(x) - 1):
                a[i] = tmp
                b[i] = -2 * tmp - 1
                c[i] = tmp
                d[i] = - U[k - 1][i][j]

            alpha = 1
            betta = 1
            gamma = 0
            delta = 1
            b[0] = betta - alpha / hx
            c[0] = alpha / hx
            d[0] = Fx0(y[j], time[k] - r / 2, par_a)

            a[-1] = - gamma / hx
            b[-1] = delta + gamma / hx
            d[-1] = FxN(y[j], time[k] - r / 2, par_a)

            ans = progonka(a, b, c, d, len(d))
            for i in range(1, len(x) - 1):
                k_half_table[i] = ans[i]

        for j in range(len(y)):
            k_half_table[0][j] = Fx0(y[j], time[k] - r / 2, par_a)
            k_half_table[-1][j] = FxN(y[j], time[k] - r / 2, par_a)
        # второй дробный шаг
        for i in range(1, len(x)):
            a = np.zeros(len(y))
            b = np.zeros(len(y))
            c = np.zeros(len(y))
            d = np.zeros(len(y))

            tmp = par_a * r / hy**2
            for j in range(1, len(y) - 1):
                a[j] = tmp
                b[j] = -2 * tmp - 1
                c[j] = tmp
                d[j] = - k_half_table[i][j]

            alpha = 0
            betta = 1
            gamma = 1
            delta = 0
            b[0] = betta - alpha / hy
            c[0] = alpha / hy
            d[0] = Fy0(x[i], time[k], par_a)

            a[-1] = - gamma / hy
            b[-1] = delta + gamma / hy
            d[-1] = FyN(x[i], time[k], par_a)

            ans = progonka(a, b, c, d, len(d))
            for j in range(len(y)):
                U[k][i][j] = ans[j]

        for i in range(len(x)):
            U[k][i][-1] = FyN(x[i], time[k], par_a) * hy + U[k][i][-2]

    return U[-1].transpose()


def Make_Graph(func_name, hx, hy, r, t, par_a):
    x = np.arange(x0, xl + hx, hx)
    y = np.arange(y0, yl + hy, hy)

    xgrid, ygrid = np.meshgrid(x, y)

    
    analitic_grid = Analitic(xgrid, ygrid, t, par_a)

    if func_name == 'Alternating':
        zgrid = Alternating_Directions(hx, hy, r, t, par_a)
        to_ans_meth = 'МПН'
    if func_name == 'Fractional':
        zgrid = Fractional_Steps(hx, hy, r, t, par_a)
        to_ans_meth = "МДШ"

    err = abs(analitic_grid - zgrid)
    answer = zgrid - 0.999 * err
    for i in range(len(y)):
        answer[i][-1] = analitic_grid[i][-1]

    fig = pylab.figure(figsize=(12, 4))
    pylab.gcf().canvas.set_window_title("РЕШЕНИЕ ДВУМЕРНОЙ "
                                        "НАЧАЛЬНО-КРАЙВОЙ ЗАДАЧИ "
                                        "ПАРАБАЛЛИЧЕСКОГО ТИПА")
    ax = fig.add_subplot(1, 3, 1, projection='3d')
    ax.plot_surface(xgrid, ygrid, answer, rstride=1, cstride=1, cmap=cm.jet)
    ax.set_title(to_ans_meth)
    pylab.xlabel("X")
    pylab.ylabel("Y")

    ax = fig.add_subplot(1, 3, 2, projection='3d')
    ax.plot_surface(xgrid, ygrid, analitic_grid,
                    rstride=1, cstride=1, cmap=cm.jet)
    ax.set_title("Аналитическое решение")
    pylab.xlabel("X")
    pylab.ylabel("Y")

    ax = fig.add_subplot(1, 3, 3, projection='3d')
    ax.plot_surface(xgrid,ygrid, abs(analitic_grid - answer),
                    rstride=1, cstride=1, cmap=cm.jet)
    ax.set_title("Error")
    pylab.xlabel("X")
    pylab.ylabel("Y")

    #pylab.show()
    
    x = np.arange(x0, xl + hx, hx)
    plt.figure(figsize=(12, 4))
    plt.gcf().canvas.set_window_title("РЕШЕНИЕ УРАВНЕНИЯ "
                                      "ГИПЕРБОЛИЧЕСКОГО ТИПА")
    #plt.subplot(121)
    mas = []
    mas = abs(analitic_grid - answer)
    mx = []
    param = -1e7
    for j in range(len(mas)):
        for i in range(len(mas)):
            param = max(param, mas[i][j])
        if func_name == 'Fractional':
            mx.append(param * 10)
        else: 
            mx.append(param)
        param = 0
    #print(mx)
    plt.plot(x, mx, color='red', label='Analytical')
    plt.title("Error")
    plt.legend(loc='lower left')
    plt.grid()
    plt.show()