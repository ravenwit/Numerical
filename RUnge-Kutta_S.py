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
from numpy import exp, log, sin, cos

'g1(y,t) function'


def g1(y1, y2, t):
    return y2


'g2(y,t) function'


def g2(y1, y2, t):
    q = 1.15
    b = 0.9
    omega0 = 2.0 / 3.0
    return -q * y2 - sin(y1) + b * cos(omega0 * t)


'Calculating c1-4 coefficient for RK'


def coef(g, y1, y2, t, tau):
    c1 = tau * g(y1, y2, t)
    c2 = tau * g(y1 + c1 / 2, y2 + c1 / 2, t + tau / 2)
    c3 = tau * g(y1 + c2 / 2, y2 + c2 / 2, t + tau / 2)
    c4 = tau * g(y1 + c3, y2 + c3, t + tau)

    return (c1 + 2 * c2 + 2 * c3 + c4) / 6


'Initialization'
t = 0  # Initial time
y1 = [0]  # at t = 0, y(0) = 0
y2 = [1]  # at t = 0, y'(0) = 1
time = [0]
tau = 0.1  # Step size

while t < 100:
    'y(i+1) = y(i) + coef(y(i), t, tau)'
    y1.append(y1[-1] + coef(g1, y1[-1], y2[-1], t, tau))
    y2.append(y2[-1] + coef(g2, y1[-1], y2[-1], t, tau))
    time.append(t)
    t += tau

plt.figure(1)
plt.plot(time, y1, 'ko')
plt.plot(time, y2, 'ro')
plt.xlabel('Time in sec.')
plt.ylabel('Magnitude')
plt.legend(['y1', 'y2'])
#
# plt.figure(2)
# plt.plot(y1, y2, 'go')
# plt.xlabel('y1')
# plt.ylabel('y2')
# plt.legend(['y2 vs y1'])
#
plt.show()