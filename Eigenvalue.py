import numpy as np
import math


def __trace(A):
    tr = 0
    n = len(A)
    for i in range(n):
        for j in range(n):
            if i == j: tr += A[i][j]
    return tr


def _get_charac_coeffs(A):
    coeffs = []
    invA = []
    n = len(A)
    I = np.array([[float(i == j) for i in range(n)] for j in range(n)])
    coeffs.append(__trace(A))
    B = A
    for i in range(1, len(A)):
        if i == len(A) - 1: invA = B - coeffs[i - 1] * I
        B = np.dot(A, B - coeffs[i - 1] * I)
        coeffs.append((1/(i+1))*__trace(B))
    invA = (1/coeffs[n-1])*invA
    return coeffs, invA


def _show_eq(Coefficient):
    n = len(Coefficient)
    equation = "x^{}".format(n)
    for i in range(n):
        equation += "{}{}*x^{}".format("+" if Coefficient[i] <= 0 else "", (-1)*Coefficient[i], n - i - 1)
    print(equation)


if __name__ == '__main__':
    line = input("\n(Write in one line separated each element by comma from left to write) \nMatrix : ")
    A = []
    for i in line.replace(" ", "").split(","):
        A.append(float(i))
    n = int(math.sqrt(len(A)))
    if n * n == len(A):
        A = np.array(A).reshape(n, n)
        coeffs, invA = _get_charac_coeffs(A)
        _show_eq(coeffs)
        print("\n", invA)
