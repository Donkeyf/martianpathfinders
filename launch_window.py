import numpy as np
from scipy.constants import gravitational_constant as G
import matplotlib.pyplot as plt
import orbitaltransformsolivia as ot

## Earth data
mass_earth = 5.9722e24 # in kg
r_earth = 6378e3 # in m

mu = G * mass_earth # in m^3/s^-2

i = 34.05256922085066
raan = 348.2164637103433
argp = 103.35486285730452
r = 600e3 + r_earth # in m
e = 0
ma = 0
T = 2 * np.pi / np.sqrt(mu) * r**(3/2)
n = 2 * np.pi / T
a = r
h = np.sqrt(mu * r)

# Create time vector for propagation
t0 = ma*T/(2*np.pi)
tf = t0 + T
t_vec = np.linspace(t0, tf, 1000)

ta, ea, ma_vec = ot.compute_anomalies(t_vec, ma, e, n)
v_r, v_n = ot.compute_orbital_velocity(h, e, ta, mu)

# Calculate perifocal coordinates by calling the function
p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta, a, e, mu, h)

# Calculate ECI coordinates by calling the function
x, y, z, dx, dy, dz = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp)

# Plot 3D trace of the orbit in Earth-centred Inertial (ECI) frame
ax = plt.figure(figsize = (8, 6)).add_subplot(projection = "3d")
ax.plot(x, y, z, color = "lightseagreen")
ax.set_title("Earth Parking Orbit in ECI Frame")
ax.view_init(azim = -145) # Rotate for a better view
ax.set(xlabel = "X (m)", ylabel = "Y (m)", zlabel = "Z (m)")
ax.grid(True)
# plt.show()

## ECEF FRAME

def get_theta(year, month, day, hours):
    # Get theta for a specific time

    # Calculate the angular velocity of Earth's rotation 
    T_E = 23*3600 + 56*60 + 4 #Earth's sidereal period in seconds (23 hrs, 56 mins, 4 seconds)
    omega_E = 2*np.pi/T_E #in rad/s

    # Calculate Julian Day Number for Epoch
    UT = hours

    J0 = 367*year - int(7*(year + int((month + 9)/12))/ 4) + int(275*month / 9) + day + 1721013.5 #in days

    JD = J0 + UT/24 #in days

    # Calculate time between Julian Day J0 and J2000
    T0 = (J0 - 2451545)/36525 #dimensionless

    # Calculate Greenwich sidereal time at 0h UT
    theta_G_0 = 100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583e-8*T0**3 #in degrees
    theta_G_0 = theta_G_0 % 360 # Wrap the angle in range 0 <= x <= 360 degrees

    # Calculate the Greenwich sidereal time
    theta_G = theta_G_0 + 360.98564724*(UT/24)
    theta_G = theta_G % 360 # Wrap the angle within 0 and 360 degrees

    # Propagate Greenwich sidereal time using Earth's rotation (vector over time)
    theta = []
    for t in t_vec:
        theta_t = np.radians(theta_G) + omega_E*(t - t0) #theta with respect to time
        theta_t = np.radians(np.degrees(theta_t) % 360)  # Wrap the angle within 0 and 360 degrees
        theta.append(theta_t)

    return theta

# Find intersection 1 - change hours until intersection is found
theta1 = get_theta(2041, 1, 31, 2.5) # UT 02:30:00

# Find intersection 2 - change hours until intersection is found, different to that of theta1
theta2 = get_theta(2041, 1, 31, 17) # UT 17:00:00

# Calculate ECEF coordinates by calling the function
x_prime1, y_prime1, z_prime1, dx_prime1, dy_prime1, dz_prime1 = ot.eci_to_ecef(x, y, z, dx, dy, dz, theta1)
x_prime2, y_prime2, z_prime2, dx_prime2, dy_prime2, dz_prime2 = ot.eci_to_ecef(x, y, z, dx, dy, dz, theta2)

## Plot 3D trace of the orbit in Earth-centred Earth-fixed (ECEF) frame
ax = plt.figure(figsize = (8, 6)).add_subplot(projection = "3d")
ax.plot(x_prime1, y_prime1, z_prime1, color = "lightseagreen")
ax.plot(x_prime2, y_prime2, z_prime2, color = "lightseagreen")
ax.set_title("Earth Parking Orbit in ECEF Frame")
ax.set(xlabel = "x' (m)", ylabel = "y' (m)", zlabel = "z' (m)")
ax.grid(True)
# plt.show()

## Groundtrace

def ra_and_dec_init(x, y, z):
    # Compute right ascension and declination
    R = np.linalg.norm([x, y, z])
    l = x / R
    m = y / R
    n2 = z / R #n2 used since n is already mean motion

    # Calculate declination in degrees
    dec_rad = np.arcsin(n2)
    dec = np.degrees(dec_rad)

    # Calculate right ascension in degrees
    if m > 0:
        ra = np.degrees(np.arccos(l / np.cos(dec_rad)))
    else:
        ra = 360 - np.degrees(np.arccos(l / np.cos(dec_rad)))

    return ra, dec

def ra_and_dec(x_prime, y_prime, z_prime):
    # Initialise empty lists for ra and dec
    ra = []
    dec = []

    for j in range(len(x)):
        alpha, delta = ra_and_dec_init(x_prime[j], y_prime[j], z_prime[j])
        ra.append(alpha)
        dec.append(delta)
    
    return ra, dec

ra1, dec1 = ra_and_dec(x_prime1, y_prime1, z_prime1)
ra2, dec2 = ra_and_dec(x_prime2, y_prime2, z_prime2)
    
# Plot groundtrace
plt.figure()
plt.scatter(ra1, dec1, s=0.5, color = "limegreen", label = "Groundtrace (UT 02:30:00)")
plt.scatter(ra2, dec2, s=0.5, color = "darkgreen", label = "Groundtrace (UT 17:00:00)")
plt.scatter(136.79, -12.39, s=40, marker="x", color = "crimson", label = "Arnhem Space Centre")
plt.xlabel("East longitude (degrees)")
plt.ylabel("Latitude (degrees)")
plt.title("Groundtrace of Earth Parking Orbit")
plt.axis("Equal")
plt.axis([0, 360, -90, 90])
plt.gca().set_aspect('equal', adjustable='box')  # For equal aspect ratio
plt.axhline(0, color='k')  # The equator
plt.grid(True)
plt.legend()
plt.show()
