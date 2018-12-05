
'''

$$$ 4th order Runge-Kutta Method for ODE $$$

&&& This is a test program to solve Chaotic 
dynamics of a pendulum. The theory is described
in the book of Computational Physics by Tao Pang
page no. 90. &&&

(c) Md. Enamul Hoque
Lecturer, Department of Physics,
SUST, Sylhet - 3114, Bangladesh.

GPU v2 License.
'''

import matplotlib.pyplot as plt
from numpy import exp, log, sin, cos, pi

'g1(y,t) function'


def g1(y1, y2, t):
    return y2


'g2(y,t) function'


def g2(y1, y2, t):
    return -(pi * pi * (y1+1) / 4)


'Calculating c1-4 coefficient for RK'


def coef(g, y1, y2, t, tau):
    c1 = tau * g(y1, y2, t)
    c2 = tau * g(y1 + c1 / 2, y2 + c1 / 2, t + tau / 2)
    c3 = tau * g(y1 + c2 / 2, y2 + c2 / 2, t + tau / 2)
    c4 = tau * g(y1 + c3, y2 + c3, t + tau)

    return (c1 + 2 * c2 + 2 * c3 + c4) / 6


def runge_kutta(ini, g1, g2):
    'Initialization'
    t = 0  # Initial time
    y1 = [ini[0]]  # at t = 0, y(0) = 0
    y2 = [ini[1]]  # at t = 0, y'(0) = 1
    time = [0]
    tau = 0.01  # Step size

    while t < 1:
        'y(i+1) = y(i) + coef(y(i), t, tau)'
        y1.append(y1[-1] + coef(g1, y1[-1], y2[-1], t, tau))
        y2.append(y2[-1] + coef(g2, y1[-1], y2[-1], t, tau))
        time.append(t)
        t += tau

    return time, y1, y2


def secant(ini, x):
    dx = 0.01
    error = 1e-6
    x1 = x + dx

    while (abs(dx) > error):
        'searching for the function value'
        time, y1x, y2 = runge_kutta([ini, x], g1, g2)
        time, y1x1, y2 = runge_kutta([ini, x1], g1, g2)

        'assigning the function value'
        fx = y1x[-1] - 1
        fx1 = y1x1[-1] - 1

        'secant method'
        d = fx1 - fx
        x2 = x1 - fx1 * (x1 - x) / d
        x = x1
        x1 = x2
        dx = x1 - x

    return x1


'searching for the adjusted initial value'
a = secant(0, 10)

'Initialization'
ini = [0, a]
print(a)
time, y1, y2 = runge_kutta(ini, g1, g2)

# plt.figure(1)
plt.plot(time, y1,)
plt.plot(time, y2,)
plt.xlabel('Time.')
plt.ylabel('y')
plt.legend(['y1', 'y2'])

# plt.figure(2)
# plt.plot(y1, y2, 'g.')
# plt.xlabel('y1')
# plt.ylabel('y2')
# plt.legend(['y2 vs y1'])

plt.show()
