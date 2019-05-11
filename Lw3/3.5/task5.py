import numpy as np
import argparse

def rectangles(interval, f, h):
    a, b = interval
    x = np.linspace(a, b, int((b - a) / h + 1))

    integral = h * np.sum([f((i + j) / 2) for i, j in zip(x, x[1:])])

    return integral


def trapezium(interval, f, h):
    a, b = interval
    x = np.linspace(a, b, int((b - a) / h + 1))
    y = [f(i) for i in x]
    integral = h * sum([i + j for i, j in zip(y[1:], y)]) / 2

    return integral


def simpson(interval, f, size, h):
    a, b = interval
    n = int((b - a) / h + 1)

    first_part = 4 * np.sum([f(a + i * h) for i in range(1, n, 2)])
    second_part = 2 * np.sum([f(a + i * h) for i in range(2, n - 1, 2)])
    integral = h * (f(a) + first_part + second_part + f(b)) / 3

    return integral


def read_interval():
    return list(map(float, input().split()))


def read_h():
    return list(map(float, input().split()))


def runge_romberg(I_h, I_2h, r, p):
    return (I_h - I_2h) / (1  - r ** p) + I_h

def main():
    wolfram_val = 0.1454
    interval = [-2, 2]
    h_series = [1, 0.5]
    r = h_series[-1] / h_series[0]
    n = 2
    func = lambda x: x**(2) / (x**(3) - 27)
    rect = []
    trap = []
    simp = []

    for h in h_series:
        rectangles_val = rectangles(interval, func, h)
        rect.append(rectangles_val)
        trapezium_val = trapezium(interval, func, h)
        trap.append(trapezium_val)
        simpson_val = simpson(interval, func, n, h)
        simp.append(simpson_val)

    with open('output1', 'w') as f1:
        f1.write(f'h1 = {h_series[0]}\th2 = {h_series[-1]}')
        f1.write('\n')
        f1.write('rectangles:')
        f1.write(str(rect))
        f1.write('\t')
        f1.write(str(runge_romberg(rect[0], rect[-1], r, n)))
        f1.write('\n')
        f1.write("trapezium:")
        f1.write(str(trap))
        f1.write('\t')
        f1.write(str(runge_romberg(trap[0], trap[-1], r, n)))
        f1.write('\n')
        f1.write("simpson:")
        f1.write(str(simp))
        f1.write('\t')
        f1.write(str(runge_romberg(simp[0], simp[-1], r, n)))
        f1.write('\n')
        #f1.write(f'analytical val = {wolfram_val}')

main()