import numpy as np


def gauss_elem(a):
    r = len(a)
    c = len(a[0])
    u = np.array([[float(i == j) for i in range(r)] for j in range(c)])

    for i in range(r - 1):
        k = i
        for j in range(i + 1, r):
            if abs(a[j][i]) > abs(a[k][i]):
                k = j
        if k != i:
            a[[i, k]] = a[[k, i]]
            u[[i, k]] = u[[k, i]]

        for j in range(i + 1, r):
            m = a[j][i] / a[i][i]
            for k in range(i + 1, r):
                a[j][k] -= m * a[i][k]
            for k in range(c):
                u[j][k] -= m * u[i][k]

    for i in range(r - 1, -1, -1):
        for j in range(i + 1, r):
            m = a[i][j]
            for k in range(c):
                u[i][k] -= m * u[j][k]
        m = 1 / a[i][i]
        for j in range(c):
            u[i][j] *= m
    # Returns the inverse matrix of input matrix obtained from gauss-jordan elimination
    return u


def __get_matrix():
    _var = int(input("How much variables in LE: "))
    line = input("\n(Write in one line separated each element by comma from left to write) \nCoefficient Matrix : ")
    A = []
    for i in line.replace(" ", "").split(","):
        A.append(float(i))
    line = input("Constant terms vector: ")
    b = []
    for i in line.replace(" ", "").split(","):
        b.append(float(i))
    return np.array(A).reshape(_var, _var), np.array(b).reshape(1, _var)


def __show_sol(solution):
    if float('NaN') in solution or float('Inf') in solution or -float('Inf') in solution:
        print("\n\tThere is no solution for this SLE.")
        return
    _vars = len(solution)
    print("Solution for the given SLE: \n")
    for i in range(_vars):
        print("\tx{} = {}".format(i + 1, solution[i][0]))


if __name__ == '__main__':
    A, b = __get_matrix()
    solution = np.dot(gauss_elem(A), b.T)
    __show_sol(solution)

