
import numpy as np
import copy


def gauss(A, b):
    m = len(A)
    n = len(A[0])
    uA = copy.deepcopy(A)
    ub = copy.deepcopy(b)
    while A[0][0] == 0:
        A[[0, 1]] = A[[1, 0]]

    for i in range(0, m):
        factor = A[i][i]
        for j in range(i, m):
            if j == i:
                uA[j] = A[j] / factor
                ub[0][j] /= factor
                continue
            for k in range(i, n):
                uA[j][k] -= (A[i][k] / factor) * A[j][i]
            ub[0][j] -= (b[0][i] / factor) * A[j][i]
        A = copy.deepcopy(uA)
        b = copy.deepcopy(ub)

    return A, b


A = np.array([[2, 1, -1], [-3, -1, 2], [-2, 1, 2]], dtype=float)
b = np.array([[8, -13, -3]], dtype=float)
C = np.array([[1, 2, -1], [3, 8, 9], [2, -11, 2]], dtype=float)
d =  np.array([[6, 10, -2]], dtype=float)
uA, uB = gauss(A, b)
uC, uD = gauss(C, d)
print(np.linalg.solve(uA, uB.T))
print("\n\n")
print(np.linalg.solve(uC, uD.T))
