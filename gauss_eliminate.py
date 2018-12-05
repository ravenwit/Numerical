# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
#
#
# def pprint(A):
#     n = len(A)
#     for i in range(0, n):
#         line = ""
#         for j in range(0, n+1):
#             line += str(A[i][j]) + "\t"
#             if j == n-1:
#                 line += "| "
#         print(line)
#     print("")
#
#
# def gauss(A):
#     n = len(A)
#
#     for i in range(0, n):
#         # Search for maximum in this column
#         maxEl = abs(A[i][i])
#         maxRow = i
#         for k in range(i+1, n):
#             if abs(A[k][i]) > maxEl:
#                 maxEl = abs(A[k][i])
#                 maxRow = k
#
#         # Swap maximum row with current row (column by column)
#         for k in range(i, n+1):
#             tmp = A[maxRow][k]
#             A[maxRow][k] = A[i][k]
#             A[i][k] = tmp
#
#         # Make all rows below this one 0 in current column
#         for k in range(i+1, n):
#             c = -A[k][i]/A[i][i]
#             for j in range(i, n+1):
#                 if i == j:
#                     A[k][j] = 0
#                 else:
#                     A[k][j] += c * A[i][j]
#
#     # Solve equation Ax=b for an upper triangular matrix A
#     x = [0 for i in range(n)]
#     for i in range(n-1, -1, -1):
#         x[i] = A[i][n]/A[i][i]
#         for k in range(i-1, -1, -1):
#             A[k][n] -= A[k][i] * x[i]
#     return x
#
#
# if __name__ == "__main__":
#     from fractions import Fraction
#     n = int(input())
#
#     A = [[0 for j in range(n+1)] for i in range(n)]
#
#     # Read input data
#     for i in range(0, n):
#         line = list(map(Fraction, input().split(" ")))
#         for j, el in enumerate(line):
#             A[i][j] = el
#     input()
#
#     line = input().split(" ")
#     lastLine = map(Fraction, line)
#     for i in range(0, n):
#         A[i][n] = lastLine[i]
#
#     # Print input
#     pprint(A)
#
#     # Calculate solution
#     x = gauss(A)
#
#     # Print result
#     line = "Result:\t"
#     for i in range(0, n):
#         line += str(x[i]) + "\t"
#     print(line)
#


# def GaussElim(A):
#     n=len(A)
#
#     for i in range(0,n):   #search for maximum in this column
#         maxE1=abs(A[i][i])
#         maxRow=i
#         for k in range(i+1,n):
#             if abs(A[k][i])>maxE1:
#                 maxE1=abs(A[k][i])
#                 maxRow=k
#                 # swap maximum row with current row (column by column)
#         for k in range(i,n+1):
#             tmp=A[maxRow][k]
#             A[maxRow][k]=A[i][k]
#             A[i][k]=tmp
#
#         for k in range(i+1,n):
#             c=-A[k][i]/A[i][i]
#             for j in range(i,n+1):
#                 if i==j:
#                     A[k][j]=0
#                 else:
#                     A[k][j]+=c*A[i][j]
#                     #solve equation Ax+b for an upper triangular matrix A
#     x=[0 for i in range(n)]
#     for i in range(n-1,-1,-1):
#         x[i]=A[i][n]/A[i][i]
#         for k in range(i-1,-1,-1):
#             A[k][n]-=A[k][i]*x[i]
#     return x
#
# def main():
#     A=[[2,7,3,6,2],
#             [3,3,4,4,6],
#             [6,9,5,3,3],
#             [6,9,5,3,3],
#             [4,2,1,7,5]]
#
#     mysum1 = GaussElim(A)
#     print('a) GaussElimination Solution = {:.1f}'.format(mysum1))
#
# if __name__ == "__main__":
#     main()

# import numpy as np
# A = np.array([[1, 2, 1], [3, 8, 1], [0, 4, 1]], dtype='float')
# b = np.array([2, 12, 2])
#
# Ab = np.hstack([A, b.reshape(-1, 1)])
#
# n = len(b)
#
# for i in range(n):
#     a = Ab[i]
#
#     for j in range(i + 1, n):
#         b = Ab[j]
#         m = a[i] / b[i]
#         Ab[j] = a - m * b
#
# for i in range(n - 1, -1, -1):
#     Ab[i] = Ab[i] / Ab[i, i]
#     a = Ab[i]
#
#     for j in range(i - 1, -1, -1):
#         b = Ab[j]
#         m = a[i] / b[i]
#         Ab[j] = a - m * b
#
# x = Ab[:, 3]
# print(x)


# import copy
# from fractions import Fraction
#
#
# def gauss(a, b):
#     a = copy.deepcopy(a)
#     b = copy.deepcopy(b)
#     n = len(a)
#     p = len(b[0])
#     det = 1
#     for i in range(n - 1):
#         k = i
#         for j in range(i + 1, n):
#             if abs(a[j][i]) > abs(a[k][i]):
#                 k = j
#         if k != i:
#             a[i], a[k] = a[k], a[i]
#             b[i], b[k] = b[k], b[i]
#             det = -det
#
#         for j in range(i + 1, n):
#             t = a[j][i] / a[i][i]
#             for k in range(i + 1, n):
#                 a[j][k] -= t * a[i][k]
#             for k in range(p):
#                 b[j][k] -= t * b[i][k]
#
#     # for i in range(n - 1, -1, -1):
#     #     for j in range(i + 1, n):
#     #         t = a[i][j]
#     #         for k in range(p):
#     #             a[i][k] -= t * a[j][k]
#     #     t = 1 / a[i][i]
#     #     det *= a[i][i]
#     #     for j in range(p):
#     #         a[i][j] *= t
#
#     for i in range(n - 1, -1, -1):
#         for j in range(i + 1, n):
#             t = a[i][j]
#             for k in range(p):
#                 b[i][k] -= t * b[j][k]
#         t = 1 / a[i][i]
#         det *= a[i][i]
#         for j in range(p):
#             b[i][j] *= t
#     return det, b
#
#
# def zeromat(p, q):
#     return [[0] * q for i in range(p)]
#
#
# def matmul(a, b):
#     n, p = len(a), len(a[0])
#     p1, q = len(b), len(b[0])
#     if p != p1:
#         raise ValueError("Incompatible dimensions")
#     c = zeromat(n, q)
#     for i in range(n):
#         for j in range(q):
#             c[i][j] = sum(a[i][k] * b[k][j] for k in range(p))
#     return c
#
#
# def mapmat(f, a):
#     return [list(map(f, v)) for v in a]
#
#
# def ratmat(a):
#     return mapmat(Fraction, a)
#
#
# # As an example, compute the determinant and inverse of 3x3 magic square
#
# a = [[1, 2, 1], [3, 0, 1], [1, 4, 0]]
# b = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
# rt=[[2], [12], [2]]
# det, c = gauss(a, b)
# print(c)
# print(matmul(c, rt))


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
