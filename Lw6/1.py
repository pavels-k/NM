# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 19:31:58 2018
@author: stesn
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
pd.set_option('precision',7)
pd.set_option('display.expand_frame_repr', False)


l = 1
t = 0.0
T = 4

def solve(A,B):
    n = A.shape[1]
    P = [-A[0,1]/A[0,0]]
    Q = [B[0]/A[0,0]]
    for i in range(1, n-1):
        P.append(-A[i,i+1]/(A[i,i] + A[i,i-1]*P[i-1]))
    for i in range(1, n):
        Q.append((B[i] - A[i,i-1]*Q[i-1])/(A[i,i] + A[i,i-1]*P[i-1]))
    P.append(0)
    X = np.zeros(shape=(n))
    X[n-1] = Q[n-1]
    for i in range(n-2,-1,-1):
        X[i] = Q[i] + P[i]*X[i+1]
    return X



def yavnaya_shema(p, s,h,tau):
    N = int(l/h)
    M = int(T/tau + 1)
    Y_0 = np.zeros(shape = (N+1))
    Y_1 = np.zeros(shape = (N+1))
    for i in range(N+1):
        x =i*h
        Y_0[i] = np.exp(2*x)
    start = 2
    if(s == 1):
        Y_1 = Y_0.copy()
    elif(s == 2):
        for i in range(N+1):
            x =i*h
            Y_1[i] = np.exp(2*x) - np.exp(2*x)*tau*tau/2
    Y = np.zeros(shape = (N+1))
    Y_00 = np.zeros(shape = (N+1))
    if(p == 3):
        for i in range(1,N):
            Y[i] = tau*tau*((Y_1[i+1] - 2*Y_1[i] + Y_1[i-1])/(h*h)+ 5*Y_1[i])\
                    + 2*Y_1[i] - Y_0[i]
        Y[0] = (4*Y[1] - Y[2])/(4*h + 3)
        Y[N] = (Y[N-2] - 4*Y[N-1])/(4*h - 3)
        Y_00 = Y_0
        Y_0 = Y_1
        Y_1 = Y.copy()
        start = 3
    for j in range(start, M+1): 
        for i in range(1,N):
            Y[i] = tau*tau*((Y_1[i+1] - 2*Y_1[i] + Y_1[i-1])/(h*h) - 5*Y_1[i])\
                    + 2*Y_1[i] - Y_0[i]
        if(p == 1):
            Y[0] = Y[1]/(2*h + 1)
            Y[N] = Y[N-1]/(1 - 2*h)
        elif(p == 2):
            Y[0] = (4*Y[1] - Y[2])/(4*h + 3)
            Y[N] = (Y[N-2] - 4*Y[N-1])/(4*h - 3)
        elif(p == 3):
            Y[0] = (Y[1] - (h*h/(2*tau*tau))*(-5*Y_1[0] + 4*Y_0[0] - Y_00[0]))/(1 + 2*h + 5*h*h/2 + h*h/(tau*tau))
            Y[N] = (Y[N-1] - (h*h/(2*tau*tau))*(-5*Y_1[N] + 4*Y_0[N] - Y_00[N]))/(1 - 2*h + 5*h*h/2 + h*h/(tau*tau))
            
        Y_00 = Y_0.copy()
        Y_0 = Y_1.copy()
        Y_1 = Y.copy()
    return Y

def neyavnaya_shema(p,s,h,tau):
    N = int(l/h)
    M = int(T/tau + 1)
    Y_0 = np.zeros(shape = (N+1))
    Y_1 = np.zeros(shape = (N+1))
    for i in range(N+1):
        x =i*h
        Y_0[i] = np.exp(2*x)
    start = 2
    if(s == 1):
        Y_1 = Y_0.copy()
    elif(s == 2):
        for i in range(N+1):
            x =i*h
            Y_1[i] = np.exp(2*x) - np.exp(2*x)*tau*tau/2
    Y = np.zeros(shape = (N+1))
    Y_00 = np.zeros(shape = (N+1))
    if(p == 3):
        for i in range(1,N):
            Y[i] = tau*tau*((Y_1[i+1] - 2*Y_1[i] + Y_1[i-1])/(h*h)+ 5*Y_1[i])\
                    + 2*Y_1[i] - Y_0[i]
        Y[0] = (4*Y[1] - Y[2])/(4*h + 3)
        Y[N] = (Y[N-2] - 4*Y[N-1])/(4*h - 3)
        Y_00 = Y_0
        Y_0 = Y_1
        Y_1 = Y.copy()
        start = 3
    
    
    for j in range(start, M+1): 
        A = np.zeros(shape = (N+1, N+1))
        B = np.zeros(shape = (N+1,))
        if(p == 1):
            A[0,0] = 2*h + 1
            A[0,1] = -1
            B[0] = 0
            A[N,N-1] = -1
            A[N,N] = 1 - 2*h
            B[N] = 0
        elif(p == 2):
            A[0,0] = 4*h + 3
            A[0,1] = -4
            A[0,2] = 1
            B[0] = 0
            A[N,N-2] = 1
            A[N,N-1] = -4
            A[N,N] = 3 - 4*h
            B[N] = 0
        elif(p == 3):
            A[0,0] = 1 + 2*h + 5*h*h/2 + h*h/(tau*tau)
            A[0,1] = -1
            B[0] = (-h*h/(2*tau*tau))*(-5*Y_1[0] + 4*Y_0[0] - Y_00[0])
            A[N,N-1] = -1
            A[N,N] = 1 - 2*h + 5*h*h/2 + h*h/(tau*tau)
            B[N] = (-h*h/(2*tau*tau))*(-5*Y_1[N] + 4*Y_0[N] - Y_00[N])
            
        for i in range(1,N):
            A[i,i-1] = -1/(h*h)
            A[i,i] = 1/(tau*tau) + 2/(h*h) + 5
            A[i,i+1] = -1/(h*h)
            B[i] = 2*Y_1[i]/(tau*tau) - Y_0[i]/(tau*tau)
            
        if(A[0,2] != 0):      
            c = (A[0,2]/A[1,2])
            A[0,:] = A[0,:] - A[1,:]*c
            B[0] = B[0] - B[1]*c
        if(A[N,N-2] != 0):
            c = (A[N,N-2]/A[N-1,N-2])
            A[N,:] = A[N,:] - A[N-1,:]*c
            B[N] = B[N] - B[N-1]*c
        
        Y = solve(A,B)
        Y_00 = Y_0.copy()
        Y_0 = Y_1.copy()
        Y_1 = Y.copy()
    return Y


# p = 1 - двухточечная первый порядок точности
# p = 2 - трехточечная второй порядок точности
# p = 3 - двухточечная второй порядок точности
#p = 3
h = l/14
tau = T/79
p = 2
s = 2
Y_1 = yavnaya_shema(p,s, h,tau)
Y_2 = neyavnaya_shema(p,s,h,tau)

U = lambda x, t : np.exp(2*x)*np.cos(t)

X = np.arange(0, l+0.01, h)
Y_true = U(X,T)
MAE_1 = np.sum(np.abs(Y_true - Y_1))/X.size
MAE_2 = np.sum(np.abs(Y_true - Y_2))/X.size
data = pd.DataFrame(data = {'x': X, 'Явная схема': Y_1, 'Неявная схема': Y_2,\
                            'Истинное':Y_true})
print('MAE явная схема:', MAE_1)
print('MAE неявная схема:', MAE_2)
print(data)
plt.figure(figsize=(10,5))
plt.plot(X,Y_1,label = 'Явная схема')
plt.plot(X,Y_2,label = 'Неявная схема')
plt.plot(X,Y_true,label = 'Аналитическое решение')
plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('u')

h_list = []
MAE_list = []
for i in range(5, 25):
    h_list.append(l/i)
    Y_3 = neyavnaya_shema(p,s,l/i,tau)
    Y_true = U(np.linspace(0, l, i + 1),T)
    MAE = np.sum(np.abs(Y_true - Y_3))/(i + 1)
    MAE_list.append(MAE)
    
plt.figure(figsize=(10,5))
plt.plot(h_list,MAE_list)
plt.xlabel('h')
plt.ylabel('MAE')
plt.title('Зависимость средней абсолютной ошибки неявной схемы от h')
tau_list = []
MAE_list = []
for i in range(30, 80):
    tau_list.append(T/i)
    Y_3 = neyavnaya_shema(p,s,h,T/i)
    Y_true = U(X,T)
    MAE = np.sum(np.abs(Y_true - Y_3))/(X.size)
    MAE_list.append(MAE)
    
plt.figure(figsize=(10,5))
plt.plot(tau_list,MAE_list)
plt.xlabel('tau')
plt.ylabel('MAE')
plt.title('Зависимость средней абсолютной ошибки неявной схемы от tau')
#hh, tt = np.meshgrid(h_list, tau_list, sparse=False, indexing='ij')
#xy = np.zeros(shape = (len(h_list), len(tau_list)))
#z_min, z_max = xy.min(), xy.max()
#for i in range(len(h_list)):
#    for j in range(len(tau_list)):
#        Y_3 = yavnaya_shema(p,s,hh[i,j],tt[i,j])
#        Y_true = U(np.linspace(0, l, int(l/hh[i,j]) + 1),T)
#        MAE = np.sum(np.abs(Y_true - Y_3))/(Y_true.size)
#        xy[i,j] = MAE
#plt.figure(figsize=(20,20))
#plt.pcolor(hh, tt, xy, vmin=0, vmax=1)
#print(np.unravel_index(xy.argmin(), xy.shape))
plt.show()