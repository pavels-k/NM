import numpy as np
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math as mth


def Analitic(x, y):
    return np.sin(x) * np.exp(y)


def Fi1(y):
    return np.exp(y)


def Fi2(y):
    return -np.exp(y)


def Fi3(x):
    return np.sin(x)


def Fi4(x):
    return mth.e * np.sin(x)


def Get_Norma(Cur, Prev):
    my_max = 0.
    for i in range(len(Cur)):
        for j in range(len(Cur[0])):
            if (abs(Cur[i][j] - Prev[i][j]) > my_max):
                my_max = abs(Cur[i][j] - Prev[i][j])
    return my_max


iterat = 0


def Simple_Iteration_Method(hx, hy, epsilon):
    x = np.arange(0, np.pi + hx, hx)
    y = np.arange(0, 1 + hy, hy)

    Cur_T = np.zeros((len(x), len(y)))

    global iterat
    iterat = 0
    for i in range(len(x)):
        Cur_T[i][0] = np.sin(x[i])
        Cur_T[i][-1] = np.sin(x[i]) * np.e
    while True:
        iterat += 1
        Prev_T = Cur_T.copy()

        for j in range(1, len(y) - 1):
            Cur_T[0][j] = Cur_T[1][j] - hx * np.exp(y[j]) # approx 1
            Cur_T[len(x) - 1][j] = Cur_T[len(x) - 2][j] - hx * np.exp(y[j])

        for i in range(1, len(Cur_T) - 1):
            for j in range(1, len(Cur_T[0]) - 1):
                Cur_T[i][j] = ((Prev_T[i + 1][j] +
                                Prev_T[i - 1][j]) / hx**2 +
                               (Prev_T[i][j + 1] +
                                Prev_T[i][j - 1]) / hy**2) *\
                    (hx**2 * hy**2 / (2 * hx**2 + 2 * hy**2))
        if Get_Norma(Cur_T, Prev_T) <= epsilon:
            break
    return Cur_T.transpose()


def Zeidel_Method(hx, hy, epsilon):
    x = np.arange(0, np.pi + hx, hx)
    y = np.arange(0, 1 + hy, hy)

    Cur_T = np.zeros((len(x), len(y)))

    global iterat
    iterat = 0
    for i in range(len(x)):
        Cur_T[i][0] = np.sin(x[i])
        Cur_T[i][-1] = np.sin(x[i]) * np.e
    while True:
        iterat += 1
        Prev_T = Cur_T.copy()

        for j in range(1, len(y) - 1):
            Cur_T[0][j] = Cur_T[1][j] - hx * np.exp(y[j]) #approx 1
            Cur_T[len(x) - 1][j] = Cur_T[len(x) - 2][j] - hx * np.exp(y[j])

        for i in range(1, len(Cur_T) - 1):
            for j in range(1, len(Cur_T[0]) - 1):
                Cur_T[i][j] = ((Prev_T[i + 1][j] +
                                Cur_T[i - 1][j]) / hx**2 +
                               (Prev_T[i][j + 1] +
                                Cur_T[i][j - 1]) / hy**2) *\
                    (hx**2 * hy**2 / (2 * hx**2 + 2 * hy**2))
        if Get_Norma(Cur_T, Prev_T) <= epsilon:
            break
    return Cur_T.transpose()


def Relaxation_Method(hx, hy, epsilon):
    x = np.arange(0, np.pi + hx, hx)
    y = np.arange(0, 1 + hy, hy)

    Cur_T = np.zeros((len(x), len(y)))

    global iterat
    iterat = 0
    # про параметр релаксации
    # 1 < relax_param
    relax_param = 1.5
    for i in range(len(x)):
        Cur_T[i][0] = np.sin(x[i])
        Cur_T[i][-1] = np.sin(x[i]) * np.e
    while True:
        iterat += 1
        Prev_T = Cur_T.copy()

        for j in range(1, len(y) - 1):
            Cur_T[0][j] = Cur_T[1][j] - hx * np.exp(y[j]) # approx 1
            Cur_T[len(x) - 1][j] = Cur_T[len(x) - 2][j] - hx * np.exp(y[j])

        for i in range(1, len(Cur_T) - 1):
            for j in range(1, len(Cur_T[0]) - 1):
                M = ((Prev_T[i + 1][j] +
                      Cur_T[i - 1][j]) / hx**2 +
                     (Prev_T[i][j + 1] +
                      Cur_T[i][j - 1]) / hy**2) *\
                    (hx**2 * hy**2 / (2 * hx**2 + 2 * hy**2))
                Cur_T[i][j] = (1 - relax_param) * \
                    Prev_T[i][j] + relax_param * M
        if Get_Norma(Cur_T, Prev_T) <= epsilon:
            break
    return Cur_T.transpose()


def Make_Graph(func_name, hx, hy, epsilon):
    x = np.arange(0, np.pi + hx, hx)
    y = np.arange(0, 1 + hy, hy)

    xgrid, ygrid = np.meshgrid(x, y)

    to_ans_meth = ""
    analitic_grid = Analitic(xgrid, ygrid)
    
    if func_name == 'simple':
        zgrid = Simple_Iteration_Method(hx, hy, epsilon)
        to_ans_meth = 'Simple Iteration Method'
    if func_name == 'zeidel':
        zgrid = Zeidel_Method(hx, hy, epsilon)
        to_ans_meth = "Zeidel's Method"
    if func_name == 'relax':
        zgrid = Relaxation_Method(hx, hy, epsilon)
        to_ans_meth = 'Relaxation Method'

    fig = pylab.figure(figsize=(12, 4))
    pylab.gcf().canvas.set_window_title("РЕШЕНИЕ УРАВНЕНИЯ "
                                        "ЭЛЛИПТИЧЕСКОГО ТИПА")
    ax = fig.add_subplot(1, 3, 1, projection='3d')
    ax.plot_surface(xgrid, ygrid, zgrid, rstride=1, cstride=1, cmap=cm.jet)
    ax.set_title(to_ans_meth + '\nit=' + str(iterat))

    ax = fig.add_subplot(1, 3, 2, projection='3d')
    ax.plot_surface(xgrid, ygrid, analitic_grid,
                    rstride=1, cstride=1, cmap=cm.jet)
    ax.set_title("Analytical solution")

    ax = fig.add_subplot(1, 3, 3, projection='3d')
    ax.plot_surface(xgrid, ygrid, abs(analitic_grid - zgrid),
                    rstride=1, cstride=1, cmap=cm.jet)
    ax.set_title("Error")

    pylab.show()
