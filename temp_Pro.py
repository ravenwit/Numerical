# import numpy as np
# from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
#
# re=[[5,6,8],[7,9,6], [8,7,7], [9,6,0],[6,0,5]]
#
# fig = plt.figure(figsize=(10,8))
# ax = fig.gca(projection='3d')
# # Legibility
# ax.set_title("Trajectories",fontsize=18)
# ax.set_xlabel("X Axis",fontsize=16)
# ax.set_ylabel("Y Axis",fontsize=16)
# ax.set_zlabel("Z Axis",fontsize=16)
# ax.text(17,15,0, "$\\uparrow\\, \\vec{B}$", color="red",fontsize=20)
# #ax.scatter(re[0,0],re[0,1],c='red') # Starting point
# ax.plot(re[:,0],re[:,1],re[:,2]) # Electron trajectory
# Axes3D.plot(re[:, 0],re[:,1],re[:,2]) # Electron trajectory
# # Final points
# #ax.scatter(re[-1,0],re[-1,1],re[-1,2],c='green', marker='>')
# fig.show()
# plt.draw()
# plt.show()

#
# """
# ============
# 3D animation
# ============
#
# A simple example of an animated plot... In 3D!
# """
# import numpy as np
# import matplotlib.pyplot as plt
# import mpl_toolkits.mplot3d.axes3d as p3
# import matplotlib.animation as animation
#
#
# def Gen_RandLine(length, dims=2):
#     """
#     Create a line using a random walk algorithm
#
#     length is the number of points for the line.
#     dims is the number of dimensions the line has.
#     """
#     lineData = np.empty((dims, length))
#     lineData[:, 0] = np.random.rand(dims)
#     for index in range(1, length):
#         # scaling the random numbers by 0.1 so
#         # movement is small compared to position.
#         # subtraction by 0.5 is to change the range to [-0.5, 0.5]
#         # to allow a line to move backwards.
#         step = ((np.random.rand(dims) - 0.5) * 0.1)
#         lineData[:, index] = lineData[:, index - 1] + step
#
#     return lineData
#
#
# def update_lines(num, dataLines, lines):
#     for line, data in zip(lines, dataLines):
#         # NOTE: there is no .set_data() for 3 dim data...
#         line.set_data(data[0:2, :num])
#         line.set_3d_properties(data[2, :num])
#     return lines
#
# # Attaching 3D axis to the figure
# fig = plt.figure()
# ax = p3.Axes3D(fig)
#
# # Fifty lines of random 3-D lines
# data = [Gen_RandLine(25, 3) for index in range(50)]
#
# # Creating fifty line objects.
# # NOTE: Can't pass empty arrays into 3d version of plot()
# lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]
#
# # Setting the axes properties
# ax.set_xlim3d([0.0, 1.0])
# ax.set_xlabel('X')
#
# ax.set_ylim3d([0.0, 1.0])
# ax.set_ylabel('Y')
#
# ax.set_zlim3d([0.0, 1.0])
# ax.set_zlabel('Z')
#
# ax.set_title('3D Test')
#
# # Creating the Animation object
# line_ani = animation.FuncAnimation(fig, update_lines, 25, fargs=(data, lines),
#                                    interval=50, blit=False)
#
# plt.show()
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import style

t = np.arange(0, 1200*0.00001, 0.00001)
x = (137**2/8)*(np.sqrt(1+(8*t/137)**2) -1)
plt.plot(t,x)
plt.show()
# fig = plt.figure()
# ax = p3.Axes3D(fig)
#
# x=[]
# y=[]
# z=[]
#
# def trajectory(i):
#     x.append(i)
#     y.append(x[-1]**2)
#     z.append(np.sin(y[-1]))
#     ax.plot(x[0:-1], y[0:-1], z[0:-1], c='blue')
#     ax.plot(1,2,3)
#     print(i)
#
# ani  = anim.FuncAnimation(fig, trajectory, 25, interval=100, blit=False)
# plt.show()
# # import RungeKutta as rk
# # rk.Butcher_tableau(3)
# # print(rk.b)
