import numpy as np
import matplotlib.pyplot as plt


import calculate_phasing_dv as phase
import nonHohmannTrajectory as nonHohmann
import orbitalTransforms as ot
import orbitalElements as oe
import FindStateVector as fsv


# Constants
mu_E = 3.986e5          # Earth gravitational parameter (km^3/s^2)
mu_M = 4.2828e4         # Mars gravitational parameter (km^3/s^2)
mu_S = 1.327e11         # Sun gravitational parameter (km^3/s^2)

r_E = 6378              # Earth radius (km)


#####  PHASING MANOUVRE DELTA V's IN PARKING  #####

# Circular parking orbit (500km altitude)
r_Ep = r_E + 1000
T_Ep = 2*np.pi/np.sqrt(mu_E) * r_Ep**(3/2)

# Time error
t_error1 = 0.01 / 100
t_error2 = 0.05 / 100
t_error3 = 0.1 / 100

# Required delta v's
dv1 = phase.phasing_deltav(T_Ep, r_Ep, mu_E, t_error1)
dv2 = phase.phasing_deltav(T_Ep, r_Ep, mu_E, t_error2)
dv3 = phase.phasing_deltav(T_Ep, r_Ep, mu_E, t_error3)
# print(dv1)
# print(dv2)
# print(dv3)







##### INTERPLANETARY TRAJECTORIES #####

# Julian day number of transfer dates
J0_E_d, T0_E_d = ot.julian_time_0(2040, 3, 23)          # depart Earth (01/07/2032)
J0_S_a, T0_S_a = ot.julian_time_0(2041, 4, 24)          # arrive Starman (23/06/2033)
J0_S_d, T0_S_d = ot.julian_time_0(2042, 7, 13)         # depart Starman (22/11/2033)
J0_M_a, T0_M_a = ot.julian_time_0(2043, 4, 17)          # arrive Mars (10/07/2034)
J0_M_d, T0_M_d = ot.julian_time_0(2043, 8, 1)         # depart Mars (08/11/2034)
J0_E_a, T0_E_a = ot.julian_time_0(2044, 7, 4)         # arrive Earth (24/10/2035)


# State vectors and orbital elements of planets on transfer dates
R_E_d, V_E_d, h_E, e_E, a_E, T_E, n_E, i_E, raan_E, argp_E, ta_E, ma_E = fsv.find_orbital_elements(J0_E_d, "Earth")
R_S_a, V_S_a, h_S, e_S, a_S, T_S, n_S, i_S, raan_S, argp_S, ta_S, ma_S = fsv.find_orbital_elements(J0_S_a, "Starman")
R_S_d, V_S_d, h_S, e_S, a_S, T_S, n_S, i_S, raan_S, argp_S, ta_S, ma_S = fsv.find_orbital_elements(J0_S_d, "Starman")
R_M_a, V_M_a, h_M, e_M, a_M, T_M, n_M, i_M, raan_M, argp_M, ta_M, ma_M = fsv.find_orbital_elements(J0_M_a, "Mars")
R_M_d, V_M_d, h_M, e_M, a_M, T_M, n_M, i_M, raan_M, argp_M, ta_M, ma_M = fsv.find_orbital_elements(J0_M_d, "Mars")
R_E_a, V_E_a, h_E, e_E, a_E, T_E, n_E, i_E, raan_E, argp_E, ta_E, ma_E = fsv.find_orbital_elements(J0_E_a, "Earth")


# Solve lambert's problem for velocities of transfer orbits at end points
## (1) from Earth to Starman rendezvous
tt_E_to_S = (J0_S_a-J0_E_d)*24*60*60
V_tf1_d, V_tf1_a = nonHohmann.lamberts_problem(R_E_d, R_S_a, tt_E_to_S, mu_S, "pro", -0.1)

## (2) from Starman to Mars
tt_S_to_M = (J0_M_a-J0_S_d)*24*60*60
V_tf2_d, V_tf2_a = nonHohmann.lamberts_problem(R_S_d, R_M_a, tt_S_to_M, mu_S, "pro", -0.1)

## (3) from Mars to Earth
tt_M_to_E = (J0_E_a-J0_M_d)*24*60*60
V_tf3_d, V_tf3_a = nonHohmann.lamberts_problem(R_M_d, R_E_a, tt_M_to_E, mu_S, "pro", -0.1)


# Calculating excess velocities for hyperbolic trajectories
v_exc_E_d = V_tf1_d - V_E_d                 # Escaping Earth
v_exc_M_a = V_tf2_a - V_M_a                 # Entering Mars
v_exc_M_d = V_tf3_d - V_M_d                 # Escaping Mars
v_exc_E_a = V_tf3_a - V_E_a                 # Entering Earth

print("depart earth vinf: " + str(np.linalg.norm(v_exc_E_d)))
print("arrive mars vinf: " + str(np.linalg.norm(v_exc_M_a)))
print("depart mars vinf: " + str(np.linalg.norm(v_exc_M_d)))
print("arrive earth vinf: " + str(np.linalg.norm(v_exc_E_a)))


# Simulating transfer trajectory
## Earth Orbit
t_vec_E = np.linspace(0, T_E, 10000)
ta_E, ea_E, ma_vec_E = ot.compute_anomalies(t_vec_E, ma_E, e_E, n_E)
p_E, q_E, w_E, dp_E, dq_E, dw_E = ot.elements_to_perifocal(ta_E, e_E, mu_S, h_E)
x_E, y_E, z_E, dx_E, dy_E, dz_E = ot.perifocal_to_eci(p_E, q_E, w_E, dp_E, dq_E, dw_E, i_E, raan_E, argp_E)

## Mars Orbit
t_vec_M = np.linspace(0, T_M, 10000)
ta_M, ea_M, ma_vec_M = ot.compute_anomalies(t_vec_M, ma_M, e_M, n_M)
p_M, q_M, w_M, dp_M, dq_M, dw_M = ot.elements_to_perifocal(ta_M, e_M, mu_S, h_M)
x_M, y_M, z_M, dx_M, dy_M, dz_M = ot.perifocal_to_eci(p_M, q_M, w_M, dp_M, dq_M, dw_M, i_M, raan_M, argp_M)

## Starman Orbit
t_vec_S = np.linspace(0, T_S, 10000)
ta_S, ea_S, ma_vec_S = ot.compute_anomalies(t_vec_S, ma_S, e_S, n_S)
p_S, q_S, w_S, dp_S, dq_S, dw_S = ot.elements_to_perifocal(ta_S, e_S, mu_S, h_S)
x_S, y_S, z_S, dx_S, dy_S, dz_S = ot.perifocal_to_eci(p_S, q_S, w_S, dp_S, dq_S, dw_S, i_S, raan_S, argp_S)
v_r_S, v_n_S = ot.compute_orbital_velocity(h_S, e_S, ta_S, mu_S)

## Transfer (1) Earth-Starman Orbit
h_tf1, e_tf1, a_tf1, T_tf1, n_tf1, i_tf1, raan_tf1, argp_tf1, ta_tf1, ma_tf1 = ot.find_orbital_elements(R_E_d, V_tf1_d, mu_S)
t_vec_tf1 = np.linspace(0, T_tf1, 10000)
ta_tf1, ea_tf1, ma_vec_tf1 = ot.compute_anomalies(t_vec_tf1, ma_tf1, e_tf1, n_tf1)
p_tf1, q_tf1, w_tf1, dp_tf1, dq_tf1, dw_tf1 = ot.elements_to_perifocal(ta_tf1, e_tf1, mu_S, h_tf1)
x_tf1, y_tf1, z_tf1, dx_tf1, dy_tf1, dz_tf1 = ot.perifocal_to_eci(p_tf1, q_tf1, w_tf1, dp_tf1, dq_tf1, dw_tf1, i_tf1, raan_tf1, argp_tf1)
v_r_tf1, v_n_tf1 = ot.compute_orbital_velocity(h_tf1, e_tf1, ta_tf1, mu_S)

## Transfer (2) Starman-Mars Orbit
h_tf2, e_tf2, a_tf2, T_tf2, n_tf2, i_tf2, raan_tf2, argp_tf2, ta_tf2, ma_tf2 = ot.find_orbital_elements(R_S_d, V_tf2_d, mu_S)
t_vec_tf2 = np.linspace(0, T_tf2, 10000)
ta_tf2, ea_tf2, ma_vec_tf2 = ot.compute_anomalies(t_vec_tf2, ma_tf2, e_tf2, n_tf2)
p_tf2, q_tf2, w_tf2, dp_tf2, dq_tf2, dw_tf2 = ot.elements_to_perifocal(ta_tf2, e_tf2, mu_S, h_tf2)
x_tf2, y_tf2, z_tf2, dx_tf2, dy_tf2, dz_tf2 = ot.perifocal_to_eci(p_tf2, q_tf2, w_tf2, dp_tf2, dq_tf2, dw_tf2, i_tf2, raan_tf2, argp_tf2)

## Transfer (3) Mars-Earth Orbit
h_tf3, e_tf3, a_tf3, T_tf3, n_tf3, i_tf3, raan_tf3, argp_tf3, ta_tf3, ma_tf3 = ot.find_orbital_elements(R_M_d, V_tf3_d, mu_S)
t_vec_tf3 = np.linspace(0, T_tf3, 10000)
ta_tf3, ea_tf3, ma_vec_tf3 = ot.compute_anomalies(t_vec_tf3, ma_tf3, e_tf3, n_tf3)
p_tf3, q_tf3, w_tf3, dp_tf3, dq_tf3, dw_tf3 = ot.elements_to_perifocal(ta_tf3, e_tf3, mu_S, h_tf3)
x_tf3, y_tf3, z_tf3, dx_tf3, dy_tf3, dz_tf3 = ot.perifocal_to_eci(p_tf3, q_tf3, w_tf3, dp_tf3, dq_tf3, dw_tf3, i_tf3, raan_tf3, argp_tf3)


# Key Indexes
## where transfer trajectories end, based on time
def arrival_time(orbit):
    if orbit == "tf1":
        dt_end = abs(t_vec_tf1 - tt_E_to_S)
        tf_end_t = np.where(dt_end == min(dt_end))[0][0]
    elif orbit == "tf2":
        dt_end = abs(t_vec_tf2 - tt_S_to_M)
        tf_end_t = np.where(dt_end == min(dt_end))[0][0]
    elif orbit == "tf3":
        dt_end = abs(t_vec_tf3 - tt_M_to_E)
        tf_end_t = np.where(dt_end == min(dt_end))[0][0]
    return tf_end_t


# Plotting
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot(x_tf1[:arrival_time("tf1")], y_tf1[:arrival_time("tf1")], z_tf1[:arrival_time("tf1")], "navy", linewidth=0.8, label="Transfer Trajectories")
ax.plot(x_tf2[:arrival_time("tf2")], y_tf2[:arrival_time("tf2")], z_tf2[:arrival_time("tf2")], "navy", linewidth=0.8)
ax.plot(x_tf3[:arrival_time("tf3")], y_tf3[:arrival_time("tf3")], z_tf3[:arrival_time("tf3")], "navy", linewidth=0.8)

ax.plot(x_S, y_S, z_S, "red", linestyle=":", linewidth=1.2, label="Starman")
ax.plot(x_M, y_M, z_M, "orangered", linestyle=":", linewidth=1.2, label="Mars")
ax.plot(x_E, y_E, z_E, "dodgerblue", linestyle=":", linewidth=1.2, label="Earth")

# ax.plot(x_tf1, y_tf1, z_tf1, "navy", linewidth=0.8)
# ax.plot(x_tf2, y_tf2, z_tf2, "navy", linewidth=0.8)
# ax.plot(x_tf3, y_tf3, z_tf3, "navy", linewidth=0.8)

# ax.plot(x_S[0:500], y_S[0:500], z_S[0:500], "navy", linewidth=0.8)
# ax.plot(x_tf1[arrival_time("tf1")], y_tf1[arrival_time("tf1")], z_tf1[arrival_time("tf1")], "bo", linewidth=0.8)
# ax.plot(x_S[0], y_S[0], z_S[0], "ro", linewidth=0.8)
# ax.plot(x_tf1[0], y_tf1[0], z_tf1[0], "bo", linewidth=0.8)


ax.legend(fontsize=8)

plt.xlabel("X (km)")
plt.ylabel("Y (km)")
ax.set_zlabel("Z (km)")

# ax.set_xlim3d(-0.7e8,0.3e8)
# ax.set_ylim3d(1.5e8,2.5e8)
# ax.set_zlim3d(-2e6,5e6)

ax.set_xlim3d(-3e8,3e8)
ax.set_ylim3d(-3e8,3e8)
ax.set_zlim3d(-3e8,3e8)


# Calculate delta v required to match Starman's orbit
## Checking relative velocity at Starman rendezvous 
v_rel_S = V_tf1_a - V_S_a          

## Plane change
dv_plane = 2*v_n_tf1[arrival_time("tf1")]*np.sin(abs(i_tf1-i_S)/2)

## Apse line rotation
vS = np.sqrt(v_r_S[0]**2 + v_n_S[0]**2)
vtf1 = np.sqrt(v_r_tf1[arrival_time("tf1")]**2 + v_n_tf1[arrival_time("tf1")]**2)

flight_angle_S = np.arctan(v_r_S[0]/v_n_S[0])
flight_angle_tf1 = np.arctan(v_r_tf1[arrival_time("tf1")]/v_n_tf1[arrival_time("tf1")])

dv_apse = np.sqrt(vS**2 + vtf1**2 - 2*vS*vtf1*np.cos(flight_angle_tf1-flight_angle_S))





##### TRAJECTORIES IN EARTH/MARS SOI #####

# Function for plotting representatively sized planets
def plot_planet(planet):
    if planet=="Mars":
        radius = 3390
        shade = 'orangered'
    elif planet=="Earth":
        radius = 6378
        shade = 'deepskyblue'
    angle1 = np.linspace(0, 2.*np.pi, 100)
    angle2 = np.linspace(0, np.pi, 100)
    planet_x = radius * np.outer(np.cos(angle1), np.sin(angle2))
    planet_y = radius * np.outer(np.sin(angle1), np.sin(angle2))
    planet_z = radius * np.outer(np.ones(np.size(angle1)), np.cos(angle2))
    ax.plot_surface(planet_x, planet_y, planet_z, color=shade, alpha=0.5)

# Function for finding periapse index based on distance
def periapse(x,y,z):
    r = np.zeros_like(x)
    for i in range(len(r)):
        r[i] = np.sqrt(x[i]**2+y[i]**2+z[i]**2)
    periapse_index = np.where(r==min(r))[0][0]
    return periapse_index


### (1) Departing Earth ###

# convert excess velocities from heliocentric frame to respective planetary frames
v_exc_E_d  = ot.helio_to_eci(v_exc_E_d, "Earth")

# find orbital elements of hyperbolic trajectory
h_E_d_hyp, e_E_d_hyp, a_E_d_hyp, i_E_d_hyp, raan_E_d_hyp, argp_E_d_hyp, ta_inf_E_d_hyp = ot.hyperbolic_orbital_elements(6378+600, v_exc_E_d, mu_E, 1)

# define circular parking orbit based on hyperbolic trajectory
v_E_d_park, x_E_d_park, y_E_d_park, z_E_d_park = ot.circular_parking_orbit(600, 0, i_E_d_hyp, raan_E_d_hyp, argp_E_d_hyp, mu_E)

# simulate hyperbolic trajectory
t_vec_E_d_hyp = np.linspace(0, 500000, 10000)
ta_E_d_hyp, ea_E_d_hype, ma_vec_E_d_hyp = ot.compute_hyperbolic_anomalies(t_vec_E_d_hyp, -np.pi*5/6, e_E_d_hyp, h_E_d_hyp, mu_E)
p_E_d_hyp, q_E_d_hyp, w_E_d_hyp, dp_E_d_hyp, dq_E_d_hyp, dw_E_d_hyp = ot.elements_to_perifocal(ta_E_d_hyp, e_E_d_hyp, mu_E, h_E_d_hyp)
x_E_d_hyp, y_E_d_hyp, z_E_d_hyp, dx_E_d_hyp, dy_E_d_hyp, dz_E_d_hyp = ot.perifocal_to_eci(p_E_d_hyp, q_E_d_hyp, w_E_d_hyp, dp_E_d_hyp, dq_E_d_hyp, dw_E_d_hyp, i_E_d_hyp, raan_E_d_hyp, argp_E_d_hyp)
vr_E_d_hyp, vn_E_d_hyp = ot.compute_orbital_velocity(h_E_d_hyp, e_E_d_hyp, ta_E_d_hyp, mu_E)

# periapse index
peri_E_d = periapse(x_E_d_hyp, y_E_d_hyp, z_E_d_hyp)

# plot trajectory
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(x_E_d_park, y_E_d_park, z_E_d_park, "darkblue", linewidth=1, linestyle=":", label="Parking Orbit")
ax.plot(x_E_d_hyp, y_E_d_hyp, z_E_d_hyp, "darkblue", linewidth=1, label="Hyperbolic Trajectory")
plot_planet("Earth")
plt.title("Departing Earth")
plt.xlabel("X (km)")
plt.ylabel("Y (km)")
ax.set_zlabel("Z (km)")
plt.legend()
ax.set_xlim3d(-3e4,3e4)
ax.set_ylim3d(-3e4,3e4)
ax.set_zlim3d(-3e4,3e4)


### (2) Arriving Mars ###

# convert excess velocities from heliocentric frame to respective planetary frames
v_exc_M_a  = ot.helio_to_eci(v_exc_M_a, "Mars")

# find orbital elements of hyperbolic trajectory
h_M_a_hyp, e_M_a_hyp, a_M_a_hyp, i_M_a_hyp, raan_M_a_hyp, argp_M_a_hyp, ta_inf_M_a_hyp = ot.hyperbolic_orbital_elements(3390+400, v_exc_M_a, mu_M, 0)

# define circular parking orbit based on hyperbolic trajectory
v_M_a_park, x_M_a_park, y_M_a_park, z_M_a_park = ot.circular_parking_orbit(400, 0, i_M_a_hyp, raan_M_a_hyp, argp_M_a_hyp, mu_M)

# simulate hyperbolic trajectory
t_vec_M_a_hyp = np.linspace(0, 500000, 10000)
ta_M_a_hyp, ea_M_a_hyp, ma_vec_M_a_hyp = ot.compute_hyperbolic_anomalies(t_vec_M_a_hyp, 0, e_M_a_hyp, h_M_a_hyp, mu_M)
p_M_a_hyp, q_M_a_hyp, w_M_a_hyp, dp_M_a_hyp, dq_M_a_hyp, dw_M_a_hyp = ot.elements_to_perifocal(ta_M_a_hyp, e_M_a_hyp, mu_M, h_M_a_hyp)
x_M_a_hyp, y_M_a_hyp, z_M_a_hyp, dx_M_a_hyp, dy_M_a_hyp, dz_M_a_hyp = ot.perifocal_to_eci(p_M_a_hyp, q_M_a_hyp, w_M_a_hyp, dp_M_a_hyp, dq_M_a_hyp, dw_M_a_hyp, i_M_a_hyp, raan_M_a_hyp, argp_M_a_hyp)
vr_M_a_hyp, vn_M_a_hyp = ot.compute_orbital_velocity(h_M_a_hyp, e_M_a_hyp, ta_M_a_hyp, mu_M)

# periapse index
peri_M_a = periapse(x_M_a_hyp, y_M_a_hyp, z_M_a_hyp)

# plot trajectory
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(x_M_a_park, y_M_a_park, z_M_a_park, "greenyellow", linewidth=1, linestyle=":", label="Initial Parking Orbit")
ax.plot(x_M_a_hyp, y_M_a_hyp, z_M_a_hyp, "greenyellow", linewidth=1, label="Arriving Hyperbolic Trajectory")
plot_planet("Mars")
plt.title("Arriving Mars")
plt.xlabel("X (km)")
plt.ylabel("Y (km)")
ax.set_zlabel("Z (km)")
plt.legend()
ax.set_xlim3d(-3e4,3e4)
ax.set_ylim3d(-3e4,3e4)
ax.set_zlim3d(-3e4,3e4)


### (3) Departing Mars ###

# convert excess velocities from heliocentric frame to respective planetary frames
v_exc_M_d  = ot.helio_to_eci(v_exc_M_d, "Mars")

# find orbital elements of hyperbolic trajectory
h_M_d_hyp, e_M_d_hyp, a_M_d_hyp, i_M_d_hyp, raan_M_d_hyp, argp_M_d_hyp, ta_inf_M_d_hyp = ot.hyperbolic_orbital_elements(3390+400, v_exc_M_d, mu_M, 1)

# define circular parking orbit based on hyperbolic trajectory
v_M_d_park, x_M_d_park, y_M_d_park, z_M_d_park = ot.circular_parking_orbit(400, 0, i_M_d_hyp, raan_M_d_hyp, argp_M_d_hyp, mu_M)

# simulate hyperbolic trajectory
t_vec_M_d_hyp = np.linspace(0, 200000, 10000)
ta_M_d_hyp, ea_M_d_hyp, ma_vec_M_d_hyp = ot.compute_hyperbolic_anomalies(t_vec_M_d_hyp, 0, e_M_d_hyp, h_M_d_hyp, mu_M)
p_M_d_hyp, q_M_d_hyp, w_M_d_hyp, dp_M_d_hyp, dq_M_d_hyp, dw_M_d_hyp = ot.elements_to_perifocal(ta_M_d_hyp, e_M_d_hyp, mu_M, h_M_d_hyp)
x_M_d_hyp, y_M_d_hyp, z_M_d_hyp, dx_M_d_hyp, dy_M_d_hyp, dz_M_d_hyp = ot.perifocal_to_eci(p_M_d_hyp, q_M_d_hyp, w_M_d_hyp, dp_M_d_hyp, dq_M_d_hyp, dw_M_d_hyp, i_M_d_hyp, raan_M_d_hyp, argp_M_d_hyp)
vr_M_d_hyp, vn_M_d_hyp = ot.compute_orbital_velocity(h_M_d_hyp, e_M_d_hyp, ta_M_d_hyp, mu_M)

# periapse index
peri_M_d = periapse(x_M_d_hyp, y_M_d_hyp, z_M_d_hyp)

# plot trajectory
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot(x_M_a_park, y_M_a_park, z_M_a_park, "greenyellow", linewidth=1, linestyle=":", label="Initial Parking Orbit")
ax.plot(x_M_a_hyp, y_M_a_hyp, z_M_a_hyp, "greenyellow", linewidth=1, label="Arriving Hyperbolic Trajectory")

ax.plot(x_M_d_park, y_M_d_park, z_M_d_park, "tomato", linewidth=1, linestyle=":", label="Parking Orbit")
ax.plot(x_M_d_hyp, y_M_d_hyp, z_M_d_hyp, "tomato", linewidth=1, label="Departing Hyperbolic Trajectory")
plot_planet("Mars")
plt.xlabel("X (km)")
plt.ylabel("Y (km)")
ax.set_zlabel("Z (km)")
plt.legend()
ax.set_xlim3d(-2e4,2e4)
ax.set_ylim3d(-2e4,2e4)
ax.set_zlim3d(-2e4,2e4)


### (4) Arriving Earth ###

# convert excess velocities from heliocentric frame to respective planetary frames
v_exc_E_a  = ot.helio_to_eci(v_exc_E_a, "Earth")

# find orbital elements of hyperbolic trajectory
h_E_a_hyp, e_E_a_hyp, a_E_a_hyp, i_E_a_hyp, raan_E_a_hyp, argp_E_a_hyp, ta_inf_E_a_hyp = ot.hyperbolic_orbital_elements(6378+600, v_exc_E_a, mu_E, 0)

# define circular parking orbit based on hyperbolic trajectory
v_E_a_park, x_E_a_park, y_E_a_park, z_E_a_park = ot.circular_parking_orbit(600, 0, i_E_a_hyp, raan_E_a_hyp, argp_E_a_hyp, mu_E)

# simulate hyperbolic trajectory
t_vec_E_a_hyp = np.linspace(0, 500000, 10000)
ta_E_a_hyp, ea_E_d_hyp, ma_vec_E_d_hyp = ot.compute_hyperbolic_anomalies(t_vec_E_a_hyp, -np.pi*5/6, e_E_a_hyp, h_E_a_hyp, mu_E)
p_E_a_hyp, q_E_a_hyp, w_E_a_hyp, dp_E_a_hyp, dq_E_a_hyp, dw_E_a_hyp = ot.elements_to_perifocal(ta_E_a_hyp, e_E_a_hyp, mu_E, h_E_a_hyp)
x_E_a_hyp, y_E_a_hyp, z_E_a_hyp, dx_E_a_hyp, dy_E_a_hyp, dz_E_a_hyp = ot.perifocal_to_eci(p_E_a_hyp, q_E_a_hyp, w_E_a_hyp, dp_E_a_hyp, dq_E_a_hyp, dw_E_a_hyp, i_E_a_hyp, raan_E_a_hyp, argp_E_a_hyp)
vr_E_a_hyp, vn_E_a_hyp = ot.compute_orbital_velocity(h_E_a_hyp, e_E_a_hyp, ta_E_a_hyp, mu_E)

# periapse index
peri_E_a = periapse(x_E_a_hyp, y_E_a_hyp, z_E_a_hyp)

# plot trajectory
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(x_E_a_hyp, y_E_a_hyp, z_E_a_hyp, "darkblue", linewidth=1, label="Arriving Earth")
plot_planet("Earth")
plt.xlabel("X (km)")
plt.ylabel("Y (km)")
ax.set_zlabel("Z (km)")
plt.legend()
ax.set_xlim3d(-3e4,3e4)
ax.set_ylim3d(-3e4,3e4)
ax.set_zlim3d(-3e4,3e4)

# plt.show()

### Delta V's Analysis ###
# assuming circular orbits everywhere

dv_E_d = vn_E_d_hyp[peri_E_d] - v_E_d_park
dv_M_a = vn_M_a_hyp[peri_M_a] - v_M_a_park
dv_M_d = vn_M_d_hyp[peri_M_d] - v_M_d_park
dv_E_a = vn_E_a_hyp[peri_E_a] - v_E_a_park

print("dv to to go from parking to hyperbolic: " + str(dv_E_d))
print("dv to to go from parking to hyperbolic: " + str(dv_M_a))
print("dv to to go from parking to hyperbolic: " + str(dv_M_d))
print("dv to to go from parking to hyperbolic: " + str(dv_E_a))

##### ERROR ANALYSIS #####
print("Mars arrival excess (MCI)"+str(v_exc_M_a))
print("Earth arrival excess (ECI)"+str(v_exc_E_a))

print(np.linalg.norm(V_tf1_a-V_S_a))
print(np.linalg.norm(V_tf2_d-V_S_d))

plt.show()