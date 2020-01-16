import numpy as np

Cur_T = np.zeros([3, 3, 3])
j = 0
k = 0
l = 0
for j in range(3):    
    for i in range(3):
        Cur_T[i][j][k] = l
        l = l + 1
Cur_T[0][2][2] = 10
Cur_T[1][2][2] = 11
Cur_T[2][2][2] = 12
print(Cur_T[-1].transpose())
