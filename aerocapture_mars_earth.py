import numpy as np
import scipy as sp
import aerocapture_master as am
import find_hyperbolic_orbital_elements as hoe
import orbitalTransforms as ot
import findOrbitalElements as foe
import matplotlib.pyplot as plt
import J2

#constants
mu_mars = 4.282831567e13
mu_earth = 3.986e14
r_earth = 6378e3
r_mars = 3396.2e3
v_inf_mars = [1.55542024e3, 3.19561118e3, 1.50728708e3]
v_inf_earth = [3.25166603, 0.07191335, 0.3165932 ]


#find the orbital elements of the hyperbolic entrance trajectory
h_m, e_m, a_m, i_m, raan_m, argp_m, ta_inf_m = hoe.hyperbolic_orbital_elements(r_mars + 24620, np.array(v_inf_mars), mu_mars, 0)

ta = np.linspace(-(ta_inf_m-0.5), ta_inf_m-0.5, 10000)
p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta, a_m, e_m, mu_mars, h_m)
x, y, z, dx, dy, dz = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i_m, raan_m, argp_m)

y0 = [x[0], y[0], z[0], dx[0], dy[0], dz[0]]

t_span = [0, 4570]
t_eval = np.arange(0, 4570, 10)

#numerically solve the aerocapture
dragpath = am.aero_capture(y0, t_span, t_eval,  27809, mu_mars, r_mars, 0)

i =  int(len(dragpath.y[0]) * 4/7)
max = 0
max_ind = 0

#find maximum and minimum altitudes
while i < len(dragpath.y[0]):
    if np.linalg.norm([dragpath.y[0][i], dragpath.y[1][i], dragpath.y[2][i]]) > max:
        max = np.linalg.norm([dragpath.y[0][i], dragpath.y[1][i], dragpath.y[2][i]])
        max_ind = i
    i += 1

min = 1000000000
j = 0
while j < len(dragpath.y[0]):
    if np.linalg.norm([dragpath.y[0][j], dragpath.y[1][j], dragpath.y[2][j]]) < min:
        min = np.linalg.norm([dragpath.y[0][j], dragpath.y[1][j], dragpath.y[2][j]])
    j += 1

# print(max - r_mars)
# print(max_ind)
# print(min - r_mars)


#############Find Parking orbit Orbital elements#################
r_a_aero = [dragpath.y[0][-1], dragpath.y[1][-1], dragpath.y[2][-1]]
v_a_aero = [dragpath.y[3][-1], dragpath.y[4][-1], dragpath.y[5][-1]]

# print(np.linalg.norm(r_a_aero))
H_1 = np.cross(r_a_aero, v_a_aero)
h_1 = np.linalg.norm(H_1)
h_unit = H_1/h_1
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
a_aero = h_1**2/mu_mars

v_aero = h_1/np.linalg.norm(r_a_aero)
v_park = h_park/np.linalg.norm(r_a_aero)
delta_v = v_park - v_aero
print(f"Delta-V Parking Orbit:{delta_v}")

ta_park = np.arange(-np.pi, np.pi, np.pi/180)
p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta_park, a_park, e_park, mu_mars, h_park)
x_m, y_m, z_m, dx_m, dy_m, dz_m = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i_park, raan_park, 0)

# print(np.linalg.norm([dx_m[0], dy_m[0], dz_m[0]]))

t = 121 * 60 * 60 * 24
raan_j2 = raan_park + J2.compute_new_argp(mu_mars,1.9555e-3 , r_mars, e_park, a_park, i_park, t)


p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta_park, a_park, e_park, mu_mars, h_park)
x_j2, y_j2, z_j2, dx_j2, dy_j2, dz_j2 = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i_park, raan_j2, 0)

r_mars = 3396.2e3
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_mars = r_mars * np.outer(np.cos(u), np.sin(v))
y_mars = r_mars * np.outer(np.sin(u), np.sin(v))
z_mars = r_mars * np.outer(np.ones(np.size(u)), np.cos(v))


############################################EARTH AEROCAPTURE####################################################################
#repeat above for earth

h_e, e_e, a_e, i_e, raan_e, argp_e, ta_inf_e = hoe.hyperbolic_orbital_elements(r_earth + 70120, np.array(v_inf_earth), mu_earth, 0)

ta_e = np.linspace(-(ta_inf_m-0.5), ta_inf_m-0.5, 10000)
p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta_e, a_e, e_e, mu_earth, h_e)
x, y, z, dx, dy, dz = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i_e, raan_e, argp_e)

y0 = [x[0], y[0], z[0], dx[0], dy[0], dz[0]]

t_span = [0, 3720]
t_eval = np.arange(0, 3720, 10)

dragpathe = am.aero_capture(y0, t_span, t_eval, 4712, mu_earth, r_earth, 1)

i =  int(len(dragpathe.y[0]) * 4/7)
max = 0
max_ind = 0

while i < len(dragpathe.y[0]):
    if np.linalg.norm([dragpathe.y[0][i], dragpathe.y[1][i], dragpathe.y[2][i]]) > max:
        max = np.linalg.norm([dragpathe.y[0][i], dragpathe.y[1][i], dragpathe.y[2][i]])
        max_ind = i
    i += 1


min = 1000000000
j = 0
while j < len(dragpathe.y[0]):
    if np.linalg.norm([dragpathe.y[0][j], dragpathe.y[1][j], dragpathe.y[2][j]]) < min:
        min = np.linalg.norm([dragpathe.y[0][j], dragpathe.y[1][j], dragpathe.y[2][j]])
    j += 1


# print(max_ind)

r_a_aero = [dragpathe.y[0][-1], dragpathe.y[1][-1], dragpathe.y[2][-1]]
v_a_aero = [dragpathe.y[3][-1], dragpathe.y[4][-1], dragpathe.y[5][-1]]

#print(np.linalg.norm(v_a_aero))
H_1 = np.cross(r_a_aero, v_a_aero)
h_1 = np.linalg.norm(H_1)
h_unit = H_1/h_1
h_park = np.sqrt(mu_earth * max)
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
a_park = h_park**2/mu_earth
a_aero = h_1**2/mu_earth

v_aero = h_1/np.linalg.norm(r_a_aero)
#print(h_park)
v_park = h_park/np.linalg.norm(r_a_aero)
#print(h_1)
delta_v = v_park - v_aero
print(f"Delta-V Parking Orbit Earth:{delta_v}")

ta_park = np.arange(-np.pi, np.pi, np.pi/180)
p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta_park, a_park, e_park, mu_earth, h_park)
x_e, y_e, z_e, dx_e, dy_e, dz_e = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i_park, raan_park, 0)

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_earth = r_earth * np.outer(np.cos(u), np.sin(v))
y_earth = r_earth * np.outer(np.sin(u), np.sin(v))
z_earth = r_earth * np.outer(np.ones(np.size(u)), np.cos(v))




###########################################PLOTTING################################

plt.figure(4)
ax = plt.axes(projection='3d')
ax.plot3D(dragpath.y[0], dragpath.y[1], dragpath.y[2], color='g', label='Aerocapture')
ax.plot3D(x_m, y_m, z_m, color='b', label='Parking Orbit')
ax.plot_surface(x_mars, y_mars, z_mars, color='orangered', alpha=0.5)
ax.set_title('Entry into Martian Sphere of Influence and Parking Orbit')
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_aspect('equal')
ax.legend(loc='upper right')
plt.show()

plt.figure(5)
ax = plt.axes(projection='3d')
ax.plot3D(dragpathe.y[0], dragpathe.y[1], dragpathe.y[2], color='m', label='Aerocapture')
ax.plot3D(x_e, y_e, z_e, color='r', label='Parking Orbit')
ax.plot_surface(x_earth, y_earth, z_earth, color='b', alpha=0.5)
ax.set_title('Earth SOI Re-entry')
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_aspect('equal')
ax.legend(loc='upper right')
plt.show()

