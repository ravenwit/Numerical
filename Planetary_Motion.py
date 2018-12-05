# --------------------------PLANETARY MOTION--------------------------#

# Gravitation is a conservative fore: E = T + V
# The total energy of the system can be expressed as two coupled 1st order ODEs:

# dr/dt = v              Where v is the velocity
# dv/dt = F/m            Where F is the force and m is the mass
# F = (G*m1*m2)/(r**2)   m1 and m2 are the mass of the sun and mars respectively

# Import necessary libraries:
import numpy as np
import matplotlib.pyplot as plt

# Radii of all planets in Astronomical Units:
rMer = 0.387  # Radius of Mercury in AU
rVen = 0.723  # Radius of Venus in AU
rEar = 1.00  # Radius of Earth in AU
rMar = 1.524  # Radius of Mars in AU
rJup = 5.203  # Radius of Jupiter in AU
rSat = 9.537  # Radius of Saturn in AU
rUra = 19.191  # Radius of Uranus in AU
rNep = 30.069  # Radius of Neptune in AU

# ----------INNER PLANETS (Mercury-->Mars)----------#

# Set parameters:
N = 687  # Mars days in a year
dt = 1.00 / N  # Time Step: Fractions of a year - 1 Mars day (i.e. 1/687)
mu = 4 * np.pi ** 2  # mu=4pi^2 is the Gravitational Parameter: mu = GM where G=6.67e-11 is the Universal Gravitational Constant and M is the mass of the body

# -----EARTH-----#

# Create an array, for all variables, of size N with all entries equal to zero:
xEar = np.zeros((N,))
yEar = np.zeros((N,))
vxEar = np.zeros((N,))
vyEar = np.zeros((N,))

# Initial Conditions:
xEar[0] = rEar  # (x0 = r, y0 = 0) in AU
vyEar[0] = np.sqrt(mu / rEar)  # (vx0 = 0, v_y0 = sqrt(mu/r)) AU/yr

# Implement Verlet Algorithm:
for k in range(0, N - 1):
    vxEar[k + 1] = vxEar[k] - (mu * xEar[k]) / (rEar ** 3) * dt
    xEar[k + 1] = xEar[k] + vxEar[k + 1] * dt
    vyEar[k + 1] = vyEar[k] - (mu * yEar[k]) / (rEar ** 3) * dt
    yEar[k + 1] = yEar[k] + vyEar[k + 1] * dt

# Plot:
E = []
c = 8000
for i in  range(len(vxEar)):

    norm_v = np.sqrt(vxEar[i]**2+vyEar[i]**2)
    beta=norm_v/c
    gama = np.sqrt(1 - beta**2)
    E_m = gama*c**2
    E_s= 3300 *c**2
    E_m = 1/2 * vxEar[i]**2 + mu/xEar[i]
    E.append(E_m)

maxE = max(E)
# E = np.array(E)/maxE
plt.plot(xEar, E)
# plt.plot(xEar, yEar, 'go')
# plt.title('Circular Orbit of Earth')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.axis('equal')
plt.show()



#
# # -----MERCURY-----#
#
# # Create an array, for all variables, of size N with all entries equal to zero:
# xMer = np.zeros((N,))
# yMer = np.zeros((N,))
# vxMer = np.zeros((N,))
# vyMer = np.zeros((N,))
#
# # Initial Conditions:
# xMer[0] = rMer  # (x0 = r, y0 = 0) in AU
# vyMer[0] = np.sqrt(mu / rMer)  # (v_x0 = 0, v_y0 = sqrt(mu/r)) AU/yr
#
# # Implement Verlet Algorithm:
# for k in range(0, N - 1):
#     vxMer[k + 1] = vxMer[k] - (mu * xMer[k]) / (rMer ** 3) * dt
#     xMer[k + 1] = xMer[k] + vxMer[k + 1] * dt
#     vyMer[k + 1] = vyMer[k] - (mu * yMer[k]) / (rMer ** 3) * dt
#     yMer[k + 1] = yMer[k] + vyMer[k + 1] * dt
#
# # -----VENUS-----#
#
# # Create an array, for all variables, of size N with all entries equal to zero:
# xVen = np.zeros((N,))
# yVen = np.zeros((N,))
# vxVen = np.zeros((N,))
# vyVen = np.zeros((N,))
#
# # Initial Conditions:
# xVen[0] = rVen  # (x0 = r, y0 = 0) in AU
# vyVen[0] = np.sqrt(mu / rVen)  # (vx0 = 0, vy0 = sqrt(mu/r)) AU/yr
#
# # Implement Verlet Algorithm:
# for k in range(0, N - 1):
#     vxVen[k + 1] = vxVen[k] - (mu * xVen[k]) / (rVen ** 3) * dt
#     xVen[k + 1] = xVen[k] + vxVen[k + 1] * dt
#     vyVen[k + 1] = vyVen[k] - (mu * yVen[k]) / (rVen ** 3) * dt
#     yVen[k + 1] = yVen[k] + vyVen[k + 1] * dt
#
# # -----MARS-----#
#
# # Create an array, for all variables, of size N with all entries equal to zero:
# xMar = np.zeros((N,))
# yMar = np.zeros((N,))
# vxMar = np.zeros((N,))
# vyMar = np.zeros((N,))
#
# # Initial Conditions:
# xMar[0] = rMar  # (x0 = r, y0 = 0) in AU
# vyMar[0] = np.sqrt(mu / rMar)  # (vx0 = 0, vy0 = sqrt(mu/r)) AU/yr
#
# # Implement Verlet Algorithm:
# for k in range(0, N - 1):
#     vxMar[k + 1] = vxMar[k] - (mu * xMar[k]) / (rMar ** 3) * dt
#     xMar[k + 1] = xMar[k] + vxMar[k + 1] * dt
#     vyMar[k + 1] = vyMar[k] - (mu * yMar[k]) / (rMar ** 3) * dt
#     yMar[k + 1] = yMar[k] + vyMar[k + 1] * dt
#
# # Plot Inner Planets:
# plt.plot(xMer, yMer, 'ro', xVen, yVen, 'yo', xEar, yEar, 'go', xMar, yMar, 'bo')
# plt.scatter(0, 0, 'r')
# plt.title('Circular Orbits of Inner Planets')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.axis('equal')
# plt.show()
#
# # ----------OUTER PLANETS (Jupiter-->Neptune)----------#
#
# # Set parameters:
# N = 59800  # Neptune days in a year
# dt = 1.00 / N  # Time Step: Fractions of a year - 1 Neptune day (i.e. 1/687)
# mu = 4 * np.pi ** 2  # mu=4pi^2 is the Gravitational Parameter: mu = GM where G=6.67e-11 is the Universal Gravitational Constant and M is the mass of the body
#
# # -----JUPITER-----#
#
# # Create an array, for all variables, of size N with all entries equal to zero:
# xJup = np.zeros((N,))
# yJup = np.zeros((N,))
# vxJup = np.zeros((N,))
# vyJup = np.zeros((N,))
#
# # Initial Conditions:
# xJup[0] = rJup  # (x0 = r, y0 = 0) in AU
# vyJup[0] = np.sqrt(mu / rJup)  # (vx0 = 0, vy0 = sqrt(mu/r)) AU/yr
#
# # Implement Verlet Algorithm:
# for k in range(0, N - 1):
#     vxJup[k + 1] = vxJup[k] - (mu * xJup[k]) / (rJup ** 3) * dt
#     xJup[k + 1] = xJup[k] + vxJup[k + 1] * dt
#     vyJup[k + 1] = vyJup[k] - (mu * yJup[k]) / (rJup ** 3) * dt
#     yJup[k + 1] = yJup[k] + vyJup[k + 1] * dt
#
# # -----SATURN-----#
#
# # Create an array, for all variables, of size N with all entries equal to zero:
# xSat = np.zeros((N,))
# ySat = np.zeros((N,))
# vxSat = np.zeros((N,))
# vySat = np.zeros((N,))
#
# # Initial Conditions:
# xSat[0] = rSat  # (x0 = r, y0 = 0) in AU
# vySat[0] = np.sqrt(mu / rSat)  # (vx0 = 0, vy0 = sqrt(mu/r)) AU/yr
#
# # Implement Verlet Algorithm:
# for k in range(0, N - 1):
#     vxSat[k + 1] = vxSat[k] - (mu * xSat[k]) / (rSat ** 3) * dt
#     xSat[k + 1] = xSat[k] + vxSat[k + 1] * dt
#     vySat[k + 1] = vySat[k] - (mu * ySat[k]) / (rSat ** 3) * dt
#     ySat[k + 1] = ySat[k] + vySat[k + 1] * dt
#
# # -----URANUS-----#
#
# # Create an array, for all variables, of size N with all entries equal to zero:
# xUra = np.zeros((N,))
# yUra = np.zeros((N,))
# vxUra = np.zeros((N,))
# vyUra = np.zeros((N,))
#
# # Initial Conditions:
# xUra[0] = rUra  # (x0 = r, y0 = 0) in AU
# vyUra[0] = np.sqrt(mu / rUra)  # (vx0 = 0, vy0 = sqrt(mu/r)) AU/yr
#
# # Implement Verlet Algorithm:
# for k in range(0, N - 1):
#     vxUra[k + 1] = vxUra[k] - (mu * xUra[k]) / (rUra ** 3) * dt
#     xUra[k + 1] = xUra[k] + vxUra[k + 1] * dt
#     vyUra[k + 1] = vyUra[k] - (mu * yUra[k]) / (rUra ** 3) * dt
#     yUra[k + 1] = yUra[k] + vyUra[k + 1] * dt
#
# # -----NEPTUNE-----#
#
# # Create an array, for all variables, of size N with all entries equal to zero:
# xNep = np.zeros((N,))
# yNep = np.zeros((N,))
# vxNep = np.zeros((N,))
# vyNep = np.zeros((N,))
#
# # Initial Conditions:
# xNep[0] = rNep  # (x0 = r, y0 = 0) in AU
# vyNep[0] = np.sqrt(mu / rNep)  # (vx0 = 0, vy0 = sqrt(mu/r)) AU/yr
#
# # Implement Verlet Algorithm:
# for k in range(0, N - 1):
#     vxNep[k + 1] = vxNep[k] - (mu * xNep[k]) / (rNep ** 3) * dt
#     xNep[k + 1] = xNep[k] + vxNep[k + 1] * dt
#     vyNep[k + 1] = vyNep[k] - (mu * yNep[k]) / (rNep ** 3) * dt
#     yNep[k + 1] = yNep[k] + vyNep[k + 1] * dt
#
# # Plot Outter Planets:
# plt.plot(xJup, yJup, 'ro', xSat, ySat, 'yo', xUra, yUra, 'bo', xNep, yNep, 'ro')
# plt.scatter(0, 0, 'r')
# plt.title('Circular Orbits of Outer Planets')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.axis('equal')
# plt.show()