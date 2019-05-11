import numpy as np
import argparse


def f(x):
    return x**(2) / (x**(3) - 27)


def simpson(a, b, h, n):
    s = f(a) + f(b)
    for i in range(1, n, 2):
        s += 4 * f(a + i * h)
    for i in range(2, n-1, 2):
        s += 2 * f(a + i * h)
    return s * h / 3


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, help='File for answer')
    args = parser.parse_args()
    a = -2
    b = 2
    h1 = 1.0
    h2 = 0.5

    x1 = np.linspace(a, b, int((b - a) / h1 + 1))
    x2 = np.linspace(a, b, int((b - a) / h2 + 1))

    y_trap1 = [f(i) for i in x1]
    y_trap2 = [f(i) for i in x2]

    rect1 = h1 * sum([f((i + j) / 2) for i, j in zip(x1, x1[1:])])
    rect2 = h2 * sum([f((i + j) / 2) for i, j in zip(x2, x2[1:])])

    trap1 = h1/2 * sum([i + j for i, j in zip(y_trap1[1:], y_trap1)])
    trap2 = h2/2 * sum([i + j for i, j in zip(y_trap2[1:], y_trap2)])

    simps1 = simpson(a, b, h1, int((b - a) / h1))
    simps2 = simpson(a, b, h2, int((b - a) / h2))
    with open(args.output, 'w') as f1:
        
         f1.write(f'Method of rectangles = {rect1}\tStep = {h1}\n')
         f1.write(f'Method of rectangles = {rect2}\tStep = {h2}\n')
         f1.write(f'Error: {round(abs(rect1 - rect2) / 3, 5)}\n')
         f1.write('='*10)
         f1.write(f'\nMethod of trapeziums = {trap1}\tStep = {h1}\n')
         f1.write(f'Method of trapeziums = {trap2}\tStep = {h2}\n')
         f1.write(f'Error: {round(abs(trap1 - trap2) / 3, 5)}\n')
         f1.write('='*10)
         f1.write(f"\nSimpson's method  = {simps1}\tStep = {h1}\n")
         f1.write(f"Simpson's method  = {simps2}\tStep = {h2}\n")
         f1.write(f'Error: {round(abs(simps1 - simps2) / 15, 5)}\n')

'''
    print(f'h = {h_series[0]}\t{r}h = {h_series[-1]}')
    print('rectangles:')
    print(rect, runge_romberg(rect[0], rect[-1], r, n))
    print("trapezium:")
    print(trap, runge_romberg(trap[0], trap[-1], r, n))
    print("simpson:")
    print(simp, runge_romberg(simp[0], simp[-1], r, n))
    print(f'analytical val = {wolfram_val}')
'''

if __name__ == "__main__":
    main()
