import numpy as np
import matplotlib.pyplot as plt


def Ux0(t):
    return 0


def Uxl(t):
    return 0


def U(x):
    return np.exp(2 * x)


def Analitic(x, t):
    return np.exp(2 * x) * np.cos(t)


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

x0 = 0
xl = 1
# t = 2
param_a = 1
param_c = -5


def autofill(x0, space_step, m, n, param_a, time_step, aprox_f=1):
    Uarray = np.zeros([n, m])

    tmp_x = x0
    for j in range(m):
        Uarray[0][j] = U(tmp_x)
        if aprox_f == 1:
            Uarray[1][j] = U(tmp_x)
        if aprox_f == 2:
            Uarray[1][j] = U(tmp_x) + \
                (param_a**2 * 4 * U(tmp_x) + param_c * U(tmp_x))\
                * time_step ** 2 / 2
        tmp_x += space_step
    return Uarray


def explicit(t, m, n, aprox, aprox_f, ans_time):
    x0 = 0
    xl = 1

    space_step = (xl - x0) / (m - 1)
    time_step = t / (n - 1)

    X = np.arange(x0, xl + space_step, space_step)

    Uarray = autofill(x0, space_step, m, n, param_a, time_step, aprox_f)

    sigma = param_a**2 * time_step**2 / space_step**2

    alpha = 1.
    betta = -2.
    gamma = 1.
    delta = -2.

    for k in range(1, n - 1):
        for j in range(1, m - 1):
            Uarray[k + 1][j] = \
                Uarray[k][j + 1] * sigma +\
                Uarray[k][j] * (-2 * sigma + 2 + param_c * (time_step**2)) + \
                Uarray[k][j - 1] * sigma - \
                Uarray[k - 1][j]
        if aprox == 1:
            Uarray[k + 1][0] = alpha * Uarray[k][1] / \
                (alpha - space_step * betta)
            Uarray[k + 1][m - 1] = gamma * Uarray[k][m - 2] / \
                (gamma + space_step * delta)
        if aprox == 2:
            Uarray[k + 1][0] = \
                (Ux0((k + 1) * time_step) +
                 alpha / 2 / space_step * Uarray[k + 1][2] -
                 2 * alpha / space_step * Uarray[k + 1][1]) /\
                (-3 * alpha / 2 / space_step + betta)
            Uarray[k + 1][m - 1] = \
                (Uxl((k + 1) * time_step) -
                 alpha / 2 / space_step * Uarray[k + 1][m - 3] +
                 2 * alpha / space_step * Uarray[k + 1][m - 2]) /\
                (3 * alpha / 2 / space_step + betta)
        if aprox == 3:
            Uarray[k + 1][0] = \
                (Ux0((k + 1) * time_step) -
                 alpha * space_step / time_step / 2 * Uarray[k][0] -
                 Uarray[k + 1][1] * alpha * 2 * param_a / space_step / 2) /\
                (alpha * (-2 * param_a / space_step / 2 -
                          space_step / time_step / 2 +
                          param_c * space_step / 2) + betta)

            Uarray[k + 1][m - 1] = \
                (Uxl((k + 1) * time_step) +
                 alpha * (space_step * Uarray[k][m - 1] / 2 / time_step +
                          2 * param_a / space_step / 2 *
                          Uarray[k + 1][m - 2])) /\
                (alpha * (2 * param_a / space_step / 2 +
                          space_step / 2 / time_step -
                          param_c * space_step / 2) + betta)
    in_array = int(ans_time / time_step)
    ans_t = in_array * time_step

    plt.figure(figsize=(12, 4))
    plt.gcf().canvas.set_window_title("РЕШЕНИЕ УРАВНЕНИЯ "
                                      "ГИПЕРБОЛИЧЕСКОГО ТИПА")
    plt.subplot(121)
    plt.plot(X, Analitic(X, ans_t), color='red', label='Analytical')
    plt.plot(X, Uarray[in_array], label='Explicit')
    plt.title("Solution")
    plt.legend(loc='lower left')
    plt.grid()

    T = np.arange(0, t + time_step, time_step)
    max_analitic_in_it_time = []

    c = 0
    for i in T:
        if i == ans_time: 
            break
        else: c = c + 1

    ans = []
    z = 0

    for k in range(len(Uarray)):
        mx = -1e7
        s = 0

        for y in range(len(Uarray[0])-1):
            #print(Uarray[k][y], Analitic(X[y], T[k]))
            mx = max(mx, abs(Uarray[k][y]-Analitic(X[y], T[k])))
            #mx = mx + ((Uarray[k][y]-Analitic(X[y], T[k]))**2)
            s = s + 1
        #mx = np.sqrt(mx)
        if k == 0: print("Ошибка", mx)
        z = z + 1
        ans.append(mx)
        mx = 0
    
    
    def f(t):
        return 0.05 * t**2*np.exp(-t**2)
    
    t = np.linspace(0,t,51)
    y = f(t)
    
    
    plt.subplot(122)
    ##plt.plot(t, y, color='red', label='Error ')
    plt.plot(T, ans, color='red', label='Error ')
    plt.title("Error")
    plt.ylabel('Error')
    plt.legend(loc='upper left')
    plt.grid()
    

    plt.show()


def implicit(t, m, n, aprox, aprox_f, ans_time):
    x0 = 0
    xl = 1

    space_step = (xl - x0) / (m - 1)
    time_step = t / (n - 1)

    X = np.arange(x0, xl + space_step, space_step)

    Uarray = autofill(x0, space_step, m, n, param_a, time_step, aprox_f)

    sigma = param_a**2 * time_step**2 / space_step**2

    alpha = 1
    betta = -2
    gamma = 1
    delta = -2

    for k in range(1, n - 1):
        a = np.zeros(m)
        b = np.zeros(m)
        c = np.zeros(m)
        d = np.zeros(m)

        for j in range(1, m - 1):
            a[j] = sigma
            b[j] = -(1 + 2 * sigma)
            c[j] = sigma
            d[j] = Uarray[k - 1][j] - \
                (param_c * time_step**2 + 2) * Uarray[k][j]
        if aprox == 1:
            b[0] = betta - alpha / space_step
            c[0] = alpha / space_step
            d[0] = Ux0((k + 1) * time_step)

            a[m - 1] = - gamma / space_step
            b[m - 1] = delta + gamma / space_step
            d[m - 1] = Uxl((k + 1) * time_step)
        elif aprox == 2:
            k0 = alpha / 2 / space_step / c[1]
            c[0] = 2 * alpha / space_step + b[1] * k0
            b[0] = (-3 * alpha / 2 / space_step + betta) + a[1] * k0
            d[0] = Ux0((k + 1) * time_step) + d[1] * k0

            k1 = -(alpha / (space_step * 2)) / a[m - 2]
            a[m - 1] = (-2 * alpha / space_step) + b[m - 2] * k1
            b[m - 1] = (3 * alpha / 2 / space_step + betta) + c[m - 2] * k1
            d[m - 1] = Uxl((k + 1) * time_step) + d[m - 2] * k1
        elif aprox == 3:
            b[0] = (alpha * (-2 * param_a / space_step / 2 -
                             space_step / time_step / 2 +
                             param_c * space_step / 2) + betta)
            c[0] = alpha * 2 * param_a / space_step / 2
            d[0] = \
                (Ux0((k + 1) * time_step) -
                 alpha * space_step / time_step / 2 * Uarray[k][0])

            a[m - 1] = -alpha * 2 * param_a / space_step / 2
            b[m - 1] = alpha * (2 * param_a / space_step / 2 +
                                space_step / time_step / 2 -
                                param_c * space_step / 2) + betta
            d[m - 1] = \
                (Uxl((k + 1) * time_step) +
                 alpha * space_step / time_step / 2 * Uarray[k][m - 1])

        Y = progonka(a, b, c, d, m)
        Uarray[k + 1] = Y

    in_array = int(ans_time / time_step)
    ans_t = in_array * time_step

    plt.figure(figsize=(12, 4))
    plt.gcf().canvas.set_window_title("РЕШЕНИЕ УРАВНЕНИЯ "
                                      "ГИПЕРБОЛИЧЕСКОГО ТИПА")
    plt.subplot(121)
    plt.plot(X, Analitic(X, ans_t), color='red', label='Analytical')
    plt.plot(X, Uarray[in_array], label='Implicit')
    plt.title("Solution")
    plt.legend(loc='lower left')
    plt.grid()
    
    T = np.arange(0, t + time_step, time_step)
    
    ans = []

    for k in range(len(Uarray)):
        mx = -1e7
        s = 0

        for y in range(len(Uarray[0])-1):
            #print(Uarray[k][y], Analitic(X[y], T[k]))
            mx = max(mx, abs(Uarray[k][y]-Analitic(X[y], T[k])))
            #mx = mx + ((Uarray[k][y]-Analitic(X[y], T[k]))**2)
            s = s + 1
        #mx = np.sqrt(mx)
        if k == 0: print("Ошибка", mx)
        ans.append(mx)
        mx = 0



    max_analitic_in_it_time = []
    def f(t):
        return ans_time * t**2*np.exp(-5*t**2) + 0.0001
    t = np.linspace(0,t,51)
    y = f(t)

    
    plt.subplot(122)
    #plt.plot(t, y, color='red', label='Error ')
    plt.plot(T, ans, color='red', label='Error ')
    plt.title("Error")
    plt.ylabel('Error')

    plt.legend(loc='upper left')
    plt.grid()
    
    plt.show()
    