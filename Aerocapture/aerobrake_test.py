import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import orbitalTransforms as ot

# #parking orbit
# mu = 4.282831567e13
# e = 0
# a = 3.7885e6
# i = 0
# raan = 0
# argp = 0
# h = np.sqrt(a * mu)
# T = np.sqrt((a**3)/(mu/(4*np.pi**2)))
# n = 2 * np.pi / T
# ma = 0

# t_sec = np.arange(0, T, 1)
# ta, ea, ma_vec = ot.compute_anomalies(t_sec, ma, e, n)

# p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta, a, e, mu, h)

# x, y, z, dx, dy, dz = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp)


# #atmosphere dip
# a_a = 3616.2e3
# e_a = 0.049776
# h_a = np.sqrt(a_a * mu)
# ma_a = 0

# T_a = np.sqrt((a_a**3)/(mu/(4*np.pi**2)))
# n_a = 2 * np.pi / T_a

# t_sec_a = np.arange(0, T_a/2, 1)
# ta_a, ea, ma_vec = ot.compute_anomalies(t_sec_a, ma_a, e_a, n_a)

# ta_air = []
# j = - 3000
# while True:
#     if 0 < ta[j] < 0.9909 or 0 > ta[j] > -0.9909:
#         ta_air.append(ta[j])
#     j += 1
#     if j == 3000:
#         break

# p_a, q_a, w_a, dp, dq, dw = ot.elements_to_perifocal(ta_air, a_a, e_a, mu, h_a)
# x_a, y_a, z_a, dx, dy, dz = ot.perifocal_to_eci(p_a, q_a, w_a, dp, dq, dw, i, raan, argp)


# #aerobrake
# y0 = [x_a[0], y_a[0], z_a[0], dx[0], dy[0], dz[0]]
# print(y0)
# def drag(t, f):
#     m = 31700
#     # x = f[0]
#     # y = f[1]
#     # z = f[2]
#     # vx = f[3]
#     # vy = f[4]
#     # vz = f[5]

#     # v = np.linalg.norm([vx, vy, vz])
#     # r = np.linalg.norm([x, y, z])

#     # ax = -mu * x/(r**3) *  - (0.0587 * v**2)/m * (vx/v)
#     # ay = -mu * y/(r**3) *  - (0.0587 * v**2)/m * (vy/v)
#     # az = -mu * z/(r**3) *  - (0.0587 * v**2)/m * (vz/v)
#     # return [vx, vy, vz, ax, ay, az]

#     r_vec = f[:3] #position
#     r_dot = f[3:] #velocity

#     r = np.linalg.norm(r_vec)
#     v = np.linalg.norm(r_dot)

#     r_ddot = -mu/(r**3) * r_vec - (0.0587 * v * r_dot)/m

#     return np.concatenate((r_dot, r_ddot))

# #dragpath = sp.integrate.RK45(drag, 0, y0, 2223, rtol=1e-6, atol=1e-6, vectorized=True)
# t_span = [0,2223] 
# t_val = np.arange(0, 2223, 1)
# dragpath = sp.integrate.solve_ivp(drag, t_span, y0, dense_output=True, rtol=1e-6, atol=1e-6)

# #Plotting Mars
# r_mars = 3396.2e3
# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 100)
# x_mars = r_mars * np.outer(np.cos(u), np.sin(v))
# y_mars = r_mars * np.outer(np.sin(u), np.sin(v))
# z_mars = r_mars * np.outer(np.ones(np.size(u)), np.cos(v))

# #parking orbit
# plt.figure(2)
# ax = plt.axes(projection='3d')
# ax.plot3D(x, y, z)
# ax.plot_surface(x_mars, y_mars, z_mars)
# ax.set_title('parking orbit around Mars 400km')
# ax.set_xlabel("x (m)")
# ax.set_ylabel("y (m)")
# ax.set_zlabel("z (m)")
# plt.show()

# plt.figure(3)
# ax = plt.axes(projection='3d')
# ax.plot3D(x_a, y_a, z_a)
# ax.plot_surface(x_mars, y_mars, z_mars)
# ax.set_title('parking orbit around Mars 400km')
# ax.set_xlabel("x (m)")
# ax.set_ylabel("y (m)")
# ax.set_zlabel("z (m)")
# plt.show()

# plt.figure(4)
# ax = plt.axes(projection='3d')
# ax.plot3D(dragpath.y[0], dragpath.y[1], dragpath.y[2])
# ax.plot_surface(x_mars, y_mars, z_mars, color='r')
# ax.set_title('')
# ax.set_xlabel("x (m)")
# ax.set_ylabel("y (m)")
# ax.set_zlabel("z (m)")
# ax.set_xlim(1e6, 3e6)
# ax.set_ylim(-3e6, -1e6)
# ax.set_aspect('equal')
# plt.show()

h = np.arange(7100, 100000, 10)
T = -23.4 - 0.00222*h
p = 0.699 * np.exp(-0.00009*h) 
rho = p/(0.1921 * (T + 273.1))

plt.figure(4)
plt.plot(h, rho)
plt.show()