import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv



x = [-5.0, -3.0, -1.0, 1.0, 3.0, 5.0]
y = [2.9442, 2.8198, 2.3562, 0.7854, 0.32175, 0.1974]


x = np.array(x)
t = np.array([[i ** j for i in x] for j in reversed(range(2))])
t_tr = np.transpose(t)
g = np.dot(t, t_tr)
second = np.dot(inv(g), np.dot(t,y))
print(second)#коэффы


x = np.array(x)
t = np.array([[i ** j for i in x] for j in reversed(range(3))])
t_tr = np.transpose(t)
g = np.dot(t, t_tr)
third = np.dot(inv(g), np.dot(t,y))
print(third)# -||-


x_vals = np.linspace(x[0], x[-1])
y_sec = [np.polyval(second, i) for i in x_vals]
y_thrd = [np.polyval(third, i) for i in x_vals]

plt.scatter(x, y, color='r')
plt.plot(x_vals, y_sec, color='b')
plt.plot(x_vals, y_thrd, color='b')
plt.show()
    
y_err = [np.polyval(second, i) for i in x]
err = sum([(y_err[idx] - i) ** 2 for idx, i in enumerate(y)])
print('error 1 = ', err)

y_err = [np.polyval(third, i) for i in x]
err = sum([(y_err[idx] - i) ** 2 for idx, i in enumerate(y)])
print('error 2 =', err)