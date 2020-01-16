import matplotlib.pyplot as plt
from numpy import *


def f(t):
    return 0.05 * t**2*exp(-5*t**2)

t = linspace(0,4,51)
y = f(t)


plt.plot(t, y, color='red', label='Error ')
plt.show()
