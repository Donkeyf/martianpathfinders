import numpy as np
from scipy.constants import gravitational_constant as G
import matplotlib.pyplot as plt

r_sun = 695700e3 # in m
mass_sun = 1988500e24 # in kg
mu_sun = G * mass_sun # in m^3/s^-2

# Orbital parameters of Earth's orbit (NASA) (NOT assuming a circular orbit)
rp_e = 147.095e9 + r_sun # in m
ra_e = 152.100e9 + r_sun # in m
e_e = 0.0167
i_e = 0.0 # Earth's inclination in ecliptic plane = 0
a_e = 149.598e9 # in m
T_e = 365.256 * 24 * 3600 # sidereal orbit period - days to seconds
raan_e = 0.0 # ??
argp_e = 0.0 # ??
ma_e = 0.0 # ??

# Calculated parameters of Earth's orbit
n_e = 2 * np.pi / T_e
h_e = np.sqrt(a_e * mu_sun * (1 - e_e**2))

# Orbital parameters of Earth's orbit (NASA) (assuming a CIRCULAR orbit)
r_ec = 149.6e9 + r_sun # in m
e_ec = 0.0
i_ec = i_e
a_ec = 149.598e9 # in m
T_ec = T_e
raan_ec = 0.0
argp_ec = 0.0
ma_ec = 0.0

# Calculated parameters of Earth's orbit (CIRCULAR)
n_ec = n_e
h_ec = np.sqrt(mu_sun * r_ec)

# Orbital parameters of Mars' orbit (NASA), assuming a circular orbit
rp_m = 206.650e9 + r_sun # perihelion in m
ra_m = 249.261e9 + r_sun # aphelion in m
e_m = 0.0935 # eccentricity
i_m = np.radians(1.848) # Mars' inclination in ecliptic plane
a_m = 227.956e9 # semi-major axis in m
T_m = 686.980 * 24 * 3600 # sidereal orbit period - days to seconds
raan_m = 0.0 # ??
argp_m = 0.0 # ??
ma_m = 0.0 # ??

# Calculated parameters of Mars' orbit
n_m = 2 * np.pi / T_m
h_m = np.sqrt(a_m * mu_sun * (1 - e_m**2))

# Orbital parameters of Starman' orbit (NASA/JPL Horizons)
JD = 2459764.0 # epoch JD
rp_s = 1.4751e11 # in m
ra_s = 2.4897e11 # in m
a_s = 1.9824e11 # in m
e_s = 0.2559073
i_s = np.radians(1.0752)
raan_s = np.radians(316.92)
argp_s = np.radians(177.75)
ma_s = 0.0 # ??
T_s = 557.165 * 24 * 3600 # in seconds

# Calculated parameters of Starman's orbit
n_s = 2 * np.pi / T_s
h_s = np.sqrt(a_s * mu_sun * (1 - e_s**2))

# Time vectors - for now - later do from epoch (t0 to tf)
t_vec_e = np.linspace(0, T_e, 1000)
t_vec_m = np.linspace(0, T_m, 1000)
t_vec_s = np.linspace(0, T_s, 1000)


################################################ Orbital Simulator ################################################

def compute_anomalies(t_vec, ma, e, n=0):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity

    # Propagate mean anomaly using mean motion (vector over time)
    ma_vec = ma + n * t_vec

    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(E):
        return E - e * np.sin(E)

    ea = np.zeros(len(ma_vec)) #initialise eccentric anomaly

    for j in range(len(ma_vec)):
        E0 = ma_vec[j] #this is the initial guess
        tol = 1

        while tol >= 1e-8:
            f = kepler_equation(E0) - ma_vec[j]
            f_prime = 1 - e * np.cos(E0)
            E1 = E0 - f/f_prime
            tol = abs(E1 - E0)
            E0 = E1
        # Calc eccentric anomaly
        ea[j] = E1

    # Calc true anomaly
    ta = 2 * np.arctan(np.tan(ea/2) * np.sqrt((1 + e)/(1 - e)))

    # Wrap anomalies to -pi:pi
    ta = (ta + np.pi) % (2 * np.pi) - np.pi
    ea = (ea + np.pi) % (2 * np.pi) - np.pi
    ma_vec = (ma_vec + np.pi) % (2 * np.pi) - np.pi
    return ta, ea, ma_vec

def compute_orbital_velocity(h, e, ta, mu):
    r = h**2/mu * 1/(1 + e * np.cos(ta))

    v_r = mu / h * e * np.sin(ta) # Radial velocity
    v_n = h/r # Normal/Tangential velocity

    return v_r, v_n

def elements_to_perifocal(ta, a, e, mu, h):
    # Calc perifocal distance in m
    r = h**2/mu * 1/(1 + e * np.cos(ta))

    # Compute perifocal coordinates
    p = r * np.cos(ta) #x-coord
    q = r * np.sin(ta) #y-coord
    w = r * 0 #z-coord is 0 because it is the orbital plane - r * 0 so it remains array of equal length

    # Compute perifocal velocities
    dp = - mu/h * np.sin(ta)
    dq = mu/h * (e + np.cos(ta))
    dw = w

    return p, q, w, dp, dq, dw

def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):
    
    # Transform coordinates to ECI frame
    x = np.zeros_like(p) #initialising
    y = np.zeros_like(p)
    z = np.zeros_like(p)

    dx = np.zeros_like(p)
    dy = np.zeros_like(p)
    dz = np.zeros_like(p)

    p_to_e = perifocal_to_eci_matrix(i, raan, argp) #calling the function defining the transformation matrix

    # Define vectors of perifocal coordinates/velocities
    perifocal_coord_vec = np.array([p, q, w])
    perifocal_vel_vec = np.array([dp, dq, dw])

    # Compute vectors of ECI coordinates/velocities
    eci_coord_vec = np.dot(p_to_e, perifocal_coord_vec) #dot product of perifocal coords with transformation matrix
    eci_vel_vec = np.dot(p_to_e, perifocal_vel_vec) #perifocal velocities to ECI velocities (via matrix multiplication)

    # Compute coordinates
    x, y, z = eci_coord_vec #ECI coords assigned to array
    dx, dy, dz = eci_vel_vec #ECI velocities assigned to array

    return x, y, z, dx, dy, dz

def perifocal_to_eci_matrix(i, raan, argp):

    # Calculate transformation matrix from perifocal to ECI frame
    c_i = np.cos(i) #rotation angles/components - capture rotation and orientation aspects of transformation
    s_i = np.sin(i) #these are computed to save space in the matrix
    c_raan = np.cos(raan)
    s_raan = np.sin(raan)
    c_argp = np.cos(argp)
    s_argp = np.sin(argp)

    p_to_e = np.array([ #3x3 matrix
        [-s_raan * c_i * s_argp + c_raan * c_argp, -s_raan * c_i * c_argp - c_raan * s_argp, s_raan * s_i],
        [c_raan * c_i * s_argp + s_raan * c_argp, c_raan * c_i * c_argp - s_raan * s_argp, -c_raan * s_i],
        [s_i * s_argp, s_i * c_argp, c_i]
    ])

    return p_to_e

ta_e, ea_e, ma_vec_e = compute_anomalies(t_vec_e, ma_e, e_e, n_e)
ta_m, ea_m, ma_vec_m = compute_anomalies(t_vec_m, ma_m, e_m, n_m)
ta_s, ea_s, ma_vec_s = compute_anomalies(t_vec_s, ma_s, e_s, n_s)

v_r_e, v_n_e = compute_orbital_velocity(h_e, e_e, ta_e, mu_sun)
v_r_m, v_n_m = compute_orbital_velocity(h_m, e_m, ta_m, mu_sun)
v_r_s, v_n_s = compute_orbital_velocity(h_s, e_s, ta_s, mu_sun)

r_e = h_e**2/mu_sun * 1/(1 + e_e * np.cos(ta_e)) #radius of Earth orbit at any point (vector)
r_m = h_e**2/mu_sun * 1/(1 + e_m * np.cos(ta_m)) #radius of Mars orbit at any point (vector)
r_s = h_e**2/mu_sun * 1/(1 + e_s * np.cos(ta_s)) #radius of Starman orbit at any point (vector)

# Calculate perifocal coordinates by calling the function
p_e, q_e, w_e, dp_e, dq_e, dw_e = elements_to_perifocal(ta_e, a_e, e_e, mu_sun, h_e)
p_m, q_m, w_m, dp_m, dq_m, dw_m = elements_to_perifocal(ta_m, a_m, e_m, mu_sun, h_m)
p_s, q_s, w_s, dp_s, dq_s, dw_s = elements_to_perifocal(ta_s, a_s, e_s, mu_sun, h_s)

# Calculate ECI coordinates by calling the function
x_e, y_e, z_e, dx_e, dy_e, dz_e = perifocal_to_eci(p_e, q_e, w_e, dp_e, dq_e, dw_e, i_e, raan_e, argp_e)
x_m, y_m, z_m, dx_m, dy_m, dz_m = perifocal_to_eci(p_m, q_m, w_m, dp_m, dq_m, dw_m, i_m, raan_m, argp_m)
x_s, y_s, z_s, dx_s, dy_s, dz_s = perifocal_to_eci(p_s, q_s, w_s, dp_s, dq_s, dw_s, i_s, raan_s, argp_s)


## Plot Earth (assume circular), Mars, and Starman (elliptical) in Helocentric ("ECI") Frame
ax = plt.figure(figsize = (8, 6)).add_subplot(projection = "3d")
ax.plot(x_e, y_e, z_e, color = "lightseagreen") # Earth orbit
ax.plot(x_m, y_m, z_m, color = "hotpink") # Earth orbit
ax.plot(x_s, y_s, z_s, color = "mediumslateblue") # Earth orbit
ax.set_title("Orbits")
# ax.view_init(azim = -145) # Rotate for a better view
ax.set(xlabel = "X (m)", ylabel = "Y (m)", zlabel = "Z (m)")
ax.grid(True)

ax.set_xlim3d(-3e11, 3e11)
ax.set_ylim3d(-3e11, 3e11)
ax.set_zlim3d(-3e11,3e11)

plt.show()
