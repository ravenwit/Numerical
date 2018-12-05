
import matplotlib.pyplot as plt
import numpy as np

a = []
b = []
c = []
F_i = []


def Butcher_tableau():
    """
        This function just inserts the value of Butcher Matrix to its
        elements a, b, c.


    :return: Returns the arrays a, b, c. Which isnt necessary. a, b and c are global
             variable.
    """

    # Referring to the global versions of a, b, c and F_i
    global a
    global b
    global c
    global F_i

    # If method equals 1 then the order is 4
    root21 = np.sqrt(21)
    c = np.array([0, 1, 1 / 2, 2 / 3, (7 - root21) / 14, (7 + root21) / 14, 1])
    a = np.array([
        [0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0],
        [3 / 8, 1 / 8, 0, 0, 0, 0, 0],
        [8 / 27, 2 / 27, 8 / 27, 0, 0, 0, 0],
        [3 * (3 * root21 - 7) / 392, 8 * (7 - root21) / 392, 48 * (7 - root21) / 392, 3 * (21 - root21) / 392, 0, 0, 0],
        [-5 * (231 + 51 * root21) / 1960, -40 * (7 + root21) / 1960, -320 * root21 / 1960,
         3 * (21 + 121 * root21) / 1960, 392 * (6 + root21) / 1960, 0, 0],
        [15 * (22 + 7 * root21) / 180, 120 / 180, 40 * (7 * root21 - 5) / 180, -63 * (3 * root21 - 2) / 180,
         -14 * (49 + 9 * root21) / 180, 70 * (7 - root21) / 180, 0]
    ])
    b = np.array([9 / 180, 0, 64 / 180, 0, 49 / 180, 49 / 180, 9 / 180])

    F_i = np.array(np.zeros(len(b)))

    return a, b, c


def __func_y1(y1, y2, t):
    """
        The y1, which is actually theta in the theory.
        The derivative of y1 returns a function y2 which is just the derivative of y1
        namely angular velocity.

        The need for y1 and t in this function is that we will later use a general
        indicator for both __func_y1 and __finc_y2. Therefore, their arguments will
        have to same.

    :param y1: y1, Theta at time t
    :param y2: y2, Angular velocity at time t
    :param t: Time
    :return: y2, Angular velocity at time t
    """
    return y2

def __func_y2(y1, y2, t):
    '''
        The y1, which is actually theta in the theory.
        The derivative of y1 returns a function y2 which is just the derivative of y1
        namely angular velocity.

        The need for y1 and t in this function is that we will later use a general
        indicator for both __func_y1 and __finc_y2. Therefore, their arguments will
        have to same.

    :param y1: y1, Theta at time t
    :param y2: y2, Angular velocity at time t
    :param t: Time
    :return: y2, Angular velocity at time t
    '''
    return (-(np.pi**2)/4)*(y1+1)


def F(func, y1, y2, t, step):
    """
        Calculating the value of F_i for a single time step.

        We are actually solving two differential equations here to obtain a solution
        for the 2nd order differential.

    :param func: Which function to use, __func_y1 or __func_y2.
    :param y1: y1, Theta at time t
    :param y2: y2, Angular velocity at time t
    :param t: Time
    :param step: h, one time step
    :return: returns nothing. But update the global variable F_i for latter use.
    """
    # Referring to the global version of F_i
    global F_i
    # Calculating according to the theory. The article has a detailed explanation.
    for i in range(len(b)):
        F_i[i] = step * func((y1 + sum([a[i][j] * F_i[j] for j in range(len(a))])),
                               (y2 + sum([a[i][j] * F_i[j] for j in range(len(a))])),
                               t + step * c[i])


def _runge_kutta(init_cond, end_boundary, step, points, further):
    Butcher_tableau()
    t = [init_cond[0]]
    y1 = [init_cond[1]]
    y2 = [init_cond[2]]
    for k in range(int((end_boundary - t[0])/step)):
        F(__func_y1, y1[-1], y2[-1], t[-1], step)
        y1.append(y1[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
        F(__func_y2, y1[-1], y2[-1], t[-1], step)
        y2.append(y2[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
        t.append(t[-1] + step)

    if further:
        index = len(t)
        for k in range(points - index):
            F(__func_y1, y1[-1], y2[-1], t[-1], step)
            y1.append(y1[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
            F(__func_y2, y1[-1], y2[-1], t[-1], step)
            y2.append(y2[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
            t.append(t[-1] + step)

    return t, y1, y2


def get_func(conditions, guess, step, points):

    zero = 1e-5
    init_cond = conditions[0:2]
    bound_cond = conditions[2:4]

    init_cond.append(guess)
    guess2 = guess + step

    while abs(guess2 - guess) >= zero:

        init_cond[2] = guess
        t , y1, y2 = _runge_kutta(init_cond, bound_cond[0], step, points, 0)
        y1_boundary_f = y1[-1] - bound_cond[1]
        init_cond[2] = guess2
        t, y1, y2 = _runge_kutta(init_cond, bound_cond[0], step, points, 0)
        y1_boundary2_f = y1[-1] - bound_cond[1]

        dnmntr = (y1_boundary2_f - y1_boundary_f)
        guess3 = guess2 - y1_boundary2_f*(guess2 - guess)/dnmntr

        guess = guess2
        guess2 = guess3

    t, y1, y2 = _runge_kutta(init_cond, bound_cond[0], step, points, 1)

    return t, y1, y2, guess3


if __name__ == '__main__':
    conditions = [0, 0, 1, 1]
    guess = 0.1
    step = 0.01
    points = 1000
    t, y1, y2, guess = get_func(conditions, guess, step, points)
    print(guess)
    plt.xlabel("Time")
    plt.ylabel("Magnitude")
    plt.plot(t, y1, label='Y1', marker='o', markersize=1)
    plt.plot(t, y2, label='Y2', marker='o', markersize=1)
    plt.plot(y1, y2, label='Y1 vs Y2', marker='o', markersize=1)
    plt.legend()
    plt.show()
