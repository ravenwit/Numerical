from warnings import warn

import numpy as np
from scipy.linalg.blas import dgemm


def linear_least_squares(a, b, residuals=False):
    """
    Return the least-squares solution to a linear matrix equation.
    Solves the equation `a x = b` by computing a vector `x` that
    minimizes the Euclidean 2-norm `|| b - a x ||^2`.  The equation may
    be under-, well-, or over- determined (i.e., the number of
    linearly independent rows of `a` can be less than, equal to, or
    greater than its number of linearly independent columns).  If `a`
    is square and of full rank, then `x` (but for round-off error) is
    the "exact" solution of the equation.
    Parameters
    ----------
    a : (M, N) array_like
        "Coefficient" matrix.
    b : (M,) array_like
        Ordinate or "dependent variable" values.
    residuals : bool
        Compute the residuals associated with the least-squares solution
    Returns
    -------
    x : (M,) ndarray
        Least-squares solution. The shape of `x` depends on the shape of
        `b`.
    residuals : int (Optional)
        Sums of residuals; squared Euclidean 2-norm for each column in
        ``b - a*x``.
    """
    if type(a) != np.ndarray or not a.flags['C_CONTIGUOUS']:
        warn('Matrix a is not a C-contiguous numpy array. The solver will create a copy, which will result' + \
             ' in increased memory usage.')

    a = np.asarray(a, order='c')
    i = dgemm(alpha=1.0, a=a.T, b=a.T, trans_b=True)
    x = np.linalg.solve(i, dgemm(alpha=1.0, a=a.T, b=b)).flatten()

    if residuals:
        return x, np.linalg.norm(np.dot(a, x) - b)
    else:
        return x


if __name__ == "__main__":
    # x = np.array([0, 1, 2, 3])
    # y = np.array([-1, 0.2, 0.9, 2.1])

    x = np.array([1, 2, 3, 4])
    y = np.array([6, 5, 7, 10])

    A = np.vstack([x, np.ones(len(x))]).T
    # A = np.vstack([np.ones(len(x)), x]).T
    A = np.asarray(A, order='c')

    m, c = linear_least_squares(A, y)

    print(m, c)

    import matplotlib.pyplot as plt

    plt.plot(x, y, 'o', label='Original data', markersize=10)
    plt.plot(x, m * x + c, 'r', label='Fitted line')
    plt.legend()
    plt.show()