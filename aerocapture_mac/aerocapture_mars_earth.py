import numpy as np
import scipy as sp
import aerocapture_master as am
import find_hyperbolic_orbital_elements as hoe
import orbitalTransforms as ot
import findOrbitalElements as foe
import matplotlib.pyplot as plt

mu_mars = 4.282831567e13
r_mars = 3396.2e3

v_inf_mars = [0.22414373, -3.47108361, 0.66968008]
v_inf_earth = [3.26109495, 0.41879044, 0.44561579]

h_m, e_m, a_m, i_m, raan_m, argp_m, ta_inf_m = hoe.hyperbolic_orbital_elements(r_mars + 189960, v_inf_mars, mu_mars, 0)

ta = np.linspace(-(ta_inf_m-1), 0, 2)
p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta, a_m, e_m, mu_mars, h_m)
x, y, z, dx, dy, dz =ot.perifocal_to_eci(p, q, w, dp, dq, dw, i_m, raan_m, argp_m)

y0 = [x[0], y[0], z[0], dx[0], dy[0], dz[0]]

t_span = [0, 8640]
t_eval = np.arange(0, 8640, 10)

dragpath = am.aero_capture(y0, t_span, t_eval, 18412, mu_mars, r_mars, 0)

i =  int(len(dragpath.y[0]) * 3/4)
max = 0
max_ind = 0

# print(np.linalg.norm([dragpath.y[0][i], dragpath.y[1][i], dragpath.y[2][i]]))

while i < len(dragpath.y[0]):
    if np.linalg.norm([dragpath.y[0][i], dragpath.y[1][i], dragpath.y[2][i]]) > max:
        max = np.linalg.norm([dragpath.y[0][i], dragpath.y[1][i], dragpath.y[2][i]])
        max_ind = i
    i += 1

# min = 1000000000
# j = 0
# while j < len(dragpath.y[0]):
#     if np.linalg.norm([dragpath.y[0][j], dragpath.y[1][j], dragpath.y[2][j]]) < min:
#         min = np.linalg.norm([dragpath.y[0][j], dragpath.y[1][j], dragpath.y[2][j]])
#     j += 1

# print(max)
# print(max- r_mars)
# print(max_ind)
# print(min - r_mars)

r_a_aero = [dragpath.y[0][-1], dragpath.y[1][-1], dragpath.y[2][-1]]
v_a_aero = [dragpath.y[3][-1], dragpath.y[4][-1], dragpath.y[5][-1]]


H = np.cross(r_a_aero, v_a_aero)
h_unit = H/np.linalg.norm(H)
h_park = np.sqrt(mu_mars * max)
H = h_unit * h_park

i_park = np.arccos(H[2]/h_park)

# Node Line
N = np.cross([0,0,1], H)
n = np.linalg.norm(N)
    
# Right Ascension of the Ascending Node
if N[1] >= 0:
    raan_park = np.arccos(N[0]/n)
elif N[1] < 0:
    raan_park = 2*np.pi - np.arccos(N[0]/n)

# Eccentricity
e_park = 0

 # Semi-major Axis
a_park = h_park**2/mu_mars

ta_park = np.arange(-np.pi, np.pi, np.pi/180)
p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta_park, a_park, e_park, mu_mars, h_park)
x_m, y_m, z_m, dx_m, dy_m, dz_m = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i_park, raan_park, 0)

r_mars = 3396.2e3
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_mars = r_mars * np.outer(np.cos(u), np.sin(v))
y_mars = r_mars * np.outer(np.sin(u), np.sin(v))
z_mars = r_mars * np.outer(np.ones(np.size(u)), np.cos(v))



plt.figure(4)
ax = plt.axes(projection='3d')
ax.plot3D(dragpath.y[0], dragpath.y[1], dragpath.y[2])
ax.plot3D(x_m, y_m, z_m, color='g')
#ax.plot3D(y0[0][10:], y0[1][10:], y0[2][10:])
ax.plot_surface(x_mars, y_mars, z_mars, color='r')
ax.set_title('')
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
# ax.set_ylim(-4e6, 4e6)
# ax.set_xlim(2e6, 4e6)
ax.set_aspect('equal')
plt.show()
