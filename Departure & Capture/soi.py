import numpy as np
import scipy as sp
import orbitalTransforms2 as ot
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import sys
import format_figures as ff
import compute_params as cp

### Jack Naylor's figure formatting setup
ff.startup_plotting(font_size=12)

########################################    CONSTANTS DEFINITIONS   ########################################
mu_e = 3.986004418 * 10**14          # Standard gravitational parameter (Earth) (m^3/s^2)
mu_s = 1.327 * 10**20              # Standard gravitational parameter (Sun) (m^3/s^2)
# mu_m = 

r_E = 6378000                      # Radius of Earth (m)
m_E = 5.974 * 10**24               # Mass of Earth (kg)
R_E = 149.6 * 10**9                 # Earth-Sun distance (m)

r_S = 696340000                    # Radius of the Sun (m)
m_S = 1.989 * 10**30               # Mass of the Sun (kg)

r_M = 3397000                      # Radius of Mars (m)
m_M = 5.974 * 10**24                # Mass of Mars (kg)
R_M = 227.9 * 10**9                 # Mars-Sun distance (m)


j2 = 1.08262668 * 10**(-3)         # J2 Perturbation Constant

########################################    OPTIONS  ########################################

do_earth_dep = 1
do_mars_arr = 1
do_mars_orbit = 0
do_mars_dep = 0
do_earth_arr = 0

########################################    DEPARTURE FROM EARTH SOI   ########################################
"""
Simulate the parking orbit and hyperbolic escape trajectory given:
    - variable departure escape velocity v_esc_earth
    - variable earth velocity v_dep_earth
    - direction of escape <= IMPORTANT (from interplanetary transfer)
    - circular parking orbit with variable parameters (i,raan,argp,a)
    - no need for J2 perturbations due to short time spent here
Output:
    - 
"""

if do_earth_dep:
    print('########################################    DEPARTURE FROM EARTH SOI   ########################################')
    # general parameters
    r_soi_E = cp.soi(R_E,m_E,m_S) # radius of Earth SOI (m)

    v_earth = np.sqrt(mu_s/R_E)
    v_tr = v_earth + 1 # CHANGE THIS - transfer orbit velocity (relative to Sun)
    v_esc_earth = v_tr - v_earth # required excess escape velocity from Earth SOI (m/s)

    alt_park_earth = 1000000 # altitude of Earth parking orbit (m)
    radius_park_earth = r_E + alt_park_earth # radius of Earth parking orbit (m)

    print(f'r_soi_E = {r_soi_E/1000.0}km')
    

    # SET parking orbit parameters HERE
    i_earth_park = cp.deg_to_rad(22)
    raan_earth_park = cp.deg_to_rad(0)
    e_earth_park = 0
    argp_earth_park = cp.deg_to_rad(0)
    ma_earth_park = cp.deg_to_rad(0)
    n_earth_park = cp.compute_n(mu_e,radius_park_earth)
    t0_earth_park = "23210.35941833"  # TODO: start epoch

    a = cp.compute_semimajoraxis(n_earth_park,mu_e)
    period = cp.compute_period(mu_e,a)
    print(f'sma={a/1000.0}km')
    print(f'period={period}')

    # Create a time series (seconds) to model the XMM orbit on (uncomment for use)
    # t = np.linspace(0,604800,10000)   # One week
    # t = np.linspace(0,7200,10000)   # Two hours
    # t = np.linspace(0,4000,10000)   # 80 minutes
    n_period = 3
    n_div = 1000
    t = np.linspace(0,n_period*period,n_div)
    dt = n_period*period/n_div

    ta, a, peri, r, v_r, v_n, r_j2 = cp.tle_orbit(t, mu_e, i_earth_park,raan_earth_park,e_earth_park,argp_earth_park,ma_earth_park,n_earth_park,dt,j2,r_E)

    # print(v_r)
    # v_x = v_r[0] * np.cos(i_earth_park)

    # Calculate the altitude
    alt = ot.magnitude(peri) - r_E

    # Plot tangential velocity
    plt.figure()
    plt.title("FairDinkum Earth Tangential velocity")
    plt.xlabel("Time (s)")
    plt.ylabel("Tangential velocity (m/s)")
    plt.plot(t,v_n,color='r')
    plt.savefig("earth_velocity.png")
    plt.show()

    # Plot altitude
    plt.figure()
    plt.title("FairDinkum Earth Altitude")
    plt.xlabel("Time (s)")
    plt.ylabel("Altitude (km)")
    plt.plot(t,alt/1000.0,color='r')
    plt.savefig("earth_altitude.png")
    plt.show()

    # Plot on 3D graph
    # Earth
    theta = np.linspace(0, 2.*np.pi, 100)
    phi = np.linspace(0, np.pi, 100)
    x = r_E/1000000.0 * np.outer(np.cos(theta), np.sin(phi))
    y = r_E/1000000.0 * np.outer(np.sin(theta), np.sin(phi))
    z = r_E/1000000.0 * np.outer(np.ones(np.size(theta)), np.cos(phi))

    r_ecef = ot.eci_to_ecef(t,r,t0_earth_park)
    print('TODO: ecef conversion including earths angular velocity causing massive change in orbit')

    # Perifocal Frame
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    cp.set_aspect_equal(ax,peri[0]/1000000.0,peri[1]/1000000.0,peri[2]/1000000.0)
    ax.plot_surface(x, y, z, color='g',alpha=0.1)
    ax.scatter(peri[0]/1000000.0,peri[1]/1000000.0,peri[2]/1000000.0,color='b',marker='.',s=0.2)
    ax.set_title("FairDinkum Earth Orbital Trace in Perifocal Frame")
    ax.set_xlabel(r'$Position, p  (\times 10^3 km)$')
    ax.set_ylabel(r'$Position, q  (\times 10^3 km)$')
    ax.set_zlabel(r'$Position, w  (\times 10^3 km)$')
    ax.ticklabel_format(style='plain')
    plt.savefig("earth_perifocal.png")
    plt.show()

    # ECI Frame
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    cp.set_aspect_equal(ax,r[0]/1000000.0,r[1]/1000000.0,r[2]/1000000.0)
    ax.plot_surface(x, y, z, color='g',alpha=0.1)
    ax.scatter(r[0]/1000000.0,r[1]/1000000.0,r[2]/1000000.0,color='b',marker='.',s=0.2)
    ax.set_title("FairDinkum Earth Orbital Trace in ECI Frame")
    ax.set_xlabel(r'$Position, x  (\times 10^3 km)$')
    ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
    ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
    ax.ticklabel_format(style='plain')
    plt.savefig("earth_eci.png")
    plt.show()

    # ECEF Frame
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    cp.set_aspect_equal(ax,r_ecef[0]/1000000.0,r_ecef[1]/1000000.0,r_ecef[2]/1000000.0)
    ax.plot_surface(x, y, z, color='g',alpha=0.1)
    ax.scatter(r_ecef[0]/1000000.0,r_ecef[1]/1000000.0,r_ecef[2]/1000000.0,c=t,marker='.',s=0.2)
    ax.set_title("FairDinkum Earth Orbital Trace in ECEF Frame")
    ax.set_xlabel(r'$Position, x  (\times 10^3 km)$')
    ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
    ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
    ax.ticklabel_format(style='plain')
    plt.savefig("earth_ecef.png")
    plt.show()


    # Convert to latitude and longitude
    long = np.arctan2(r_ecef[1],r_ecef[0])/np.pi * 180
    lat = np.arctan2(r_ecef[2],np.sqrt(r_ecef[0]**2+r_ecef[1]**2))/np.pi * 180

    plt.figure()
    plt.ylim(-90,90)
    plt.xlim(-179,180)
    plt.scatter(long,lat,color='b',marker='.',s=0.5)
    plt.title("FairDinkum Earth Ground Track")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig("earth_ground_track.png")
    plt.show()

    # J2 in ECI Frame
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    cp.set_aspect_equal(ax,r[0]/1000000.0,r[1]/1000000.0,r[2]/1000000.0)
    ax.plot_surface(x, y, z, color='g',alpha=0.1)
    ax.scatter(r[0]/1000000.0,r[1]/1000000.0,r[2]/1000000.0,color='b',marker='.',s=0.2)
    ax.set_title("FairDinkum Earth Orbital Trace in ECI Frame with J2")
    ax.set_xlabel(r'$Position, x  (\times 10^3 km)$')
    ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
    ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
    ax.ticklabel_format(style='plain')
    plt.savefig("earth_j2_eci.png")
    plt.show()

    # TODO: Change to escape plane orbit? OR assume launching directly into a parking orbit on the ecliptic plane

    ### ESCAPE

    ### Given this excess velocity needed, calculate the hyperbolic escape from Earth's sphere of influence
    r_p = radius_park_earth
    e_hyp = 1 + r_p*v_esc_earth**2/mu_e # eccentricity of hyperbola
    h_hyp = r_p * np.sqrt(v_esc_earth**2 + (2*mu_e)/r_p) # angular momentum of hyperbolic trajectory
    v_p = h_hyp/r_p # tangential velocity at periapsis of hyperbolic trajectory

    v_c = np.sqrt(mu_e/r_p) # tangential velocity at circular parking orbit
    print(f'radius_park_earth = {r_p/1000.0} km')
    print(f'v_park_earth = {v_c/1000.0} km/s')

    dv_earth_esc = v_p - v_c # delta v needed to enter hyperbolic escape trajectory
    print(f'delta v to enter escape trajectory = {dv_earth_esc}')

    # Compute hyperbola parameters
    a_hyp = (h_hyp**2/mu_e) * (1/(e_hyp**2-1))
    b_hyp = a_hyp * np.sqrt(e_hyp**2-1)
    beta_hyp = np.arccos(1/e_hyp)

    # Hyperbolic escape trajectory in 3D
    theta_hyp = np.linspace(0,(11/17)*np.pi,1000)
    p_hyp = -a_hyp * (e_hyp+np.cos(theta_hyp))/(1+e_hyp*np.cos(theta_hyp)) + r_p + a_hyp
    q_hyp = b_hyp * (np.sqrt(e_hyp**2-1)*np.sin(theta_hyp))/(1+e_hyp*np.cos(theta_hyp))
    w_hyp = np.zeros(np.shape(q_hyp))


    peri_hyp = np.row_stack((p_hyp,q_hyp))
    peri_hyp = np.transpose(peri_hyp)
    # Rotate such that escape direction parallel to Earth's heliocentric velocity vector (x-axis)
    
    """
    For non-Hohmann transfer manoeuvre: hyperbolic escape direction is not parallel to earth velocity vector (ECI x-direction)
        - has inclination (CANNOT assume launch directly into parking orbit to escape from as above)
        - has rotation (implemented below)
    """
    rot_esc_earth_deg = 0
    rot_esc_earth = cp.deg_to_rad(rot_esc_earth_deg) # positive => inwards towards sun, negative => outwards away from sun
    rot_matrix = np.array([[np.cos(np.pi+beta_hyp+rot_esc_earth),-np.sin(np.pi+beta_hyp+rot_esc_earth)],[np.sin(np.pi+beta_hyp+rot_esc_earth),np.cos(np.pi+beta_hyp+rot_esc_earth)]])
    # rot_matrix = np.array([[np.cos(np.pi),-np.sin(np.pi)],[np.sin(np.pi),np.cos(np.pi)]])
    for i in range(len(peri_hyp)):
        peri_hyp[i] = np.dot(rot_matrix,peri_hyp[i])
    peri_hyp = np.transpose(peri_hyp)
    
    # recombine with w_hyp
    peri_hyp = np.row_stack((peri_hyp[0],peri_hyp[1],w_hyp))
    
    transform = ot.perifocal_to_eci_matrix(i_earth_park,raan_earth_park,argp_earth_park)    
    # Apply transform to create ECI frame coordinates array
    r_hyp = np.dot(transform,peri_hyp)

    # ECI frame + hyperbolic escape
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    cp.set_aspect_equal(ax,r[0]/1000000.0,r[1]/1000000.0,r[2]/1000000.0)
    # cp.set_aspect_equal(ax,r_soi_E/1000000.0,r_soi_E/1000000.0,r_soi_E/1000000.0)
    ax.plot_surface(x, y, z, color='g',alpha=0.1)
    ax.scatter(r[0]/1000000.0,r[1]/1000000.0,r[2]/1000000.0,color='b',marker='.',s=0.2,label='parking orbit')
    ax.scatter(r_hyp[0]/1000000.0,r_hyp[1]/1000000.0,r_hyp[2]/1000000.0,color='r',marker='.',s=0.2,label='escape eci')
    ax.set_title(f'FairDinkum Earth Orbital Trace and Escape in ECI Frame at {rot_esc_earth_deg}deg')
    ax.set_xlabel(r'$Position, x  (\times 10^3 km)$ (direction of v_earth in Sun frame)')
    ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
    ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
    ax.ticklabel_format(style='plain')
    ax.legend()
    plt.savefig("earth_escape_eci.png")
    plt.show()


########################################   ARRIVAL TO MARS SOI   ########################################
"""
Simulate the hyperbolic capture trajectory given:
    - variable arrival escape velocity v_esc_mars
    - variable mars velocity v_arr_mars
    - circular parking orbit with variable parameters (i,raan,argp) and altitude 400m
Output:
    - 
"""
if do_mars_arr:
    print('########################################   ARRIVAL TO MARS SOI   ########################################')
    # parameters
    r_soi_M = cp.soi(R_M,m_M,m_S)  # TODO: make sure this mars-sun distance is correct for the time

    v_mars = np.sqrt(mu_s/R_M) # need to change this to 
    v_cap_mars = v_mars + 1 # CHANGE THIS - capture velocity (relative to Sun)
    v_arr_mars= v_cap_mars - v_mars # required excess velocity for Mars capture (m/s)

    print(f'r_soi_M = {r_soi_M/1000.0}km')

    




########################################   ORBIT IN MARS SOI   ########################################
"""
Simulate the Mars parking orbit given:
    - circular parking orbit with variable parameters (i,raan,argp) and altitude 400m
    - propagate J2 perturbations over 3-4 months
Output:
    - 
"""
if do_mars_orbit:
    print('########################################   ORBIT IN MARS SOI   ########################################')
    # parameters
    r_soi_M = cp.soi(R_M,m_M,m_S)
    v_esc_mars = 6000 # m/s
    print(f'r_soi_M = {r_soi_M/1000.0}km')

########################################   DEPARTURE FROM MARS SOI   ########################################
"""
Simulate the hyperbolic escape trajectory given:
    - circular parking orbit with variable parameters (i,raan,argp) and altitude 400m AFTER PROPAGATING J2
    - direction of escape
    - excess escape velocity
Output:
    - departure orbit
    - hyperbolic departure
"""
if do_mars_dep:
    print('########################################   DEPARTURE FROM MARS SOI   ########################################')
    # parameters
    r_soi_M = cp.soi(R_M,m_M,m_S) # m
    v_esc_mars = 6000 # m/s
    print(f'r_soi_M = {r_soi_M/1000.0}km')

########################################   ARRIVAL TO EARTH SOI   ########################################
"""
Simulate the hyperbolic capture trajectory given:
    - circular parking orbit with variable parameters (i,raan,argp) and altitude 400m AFTER PROPAGATING J2
    - direction of escape
    - excess escape velocity
Output:
    - departure orbit
    - hyperbolic departure
"""
if do_mars_dep:
    print('########################################   ARRIVAL TO EARTH SOI   ########################################')
    # parameters
    r_soi_M = cp.soi(R_M,m_M,m_S) # m
    v_esc_mars = 6000 # m/s
    print(f'r_soi_M = {r_soi_M/1000.0}km')


exit(0)

########################################    QUESTION 1   ########################################

# Parse Keplerian orbital parameters from a single TLE (Two-Line Element) for the XMM
with open(sys.argv[1]) as f:
    lines = f.readlines()

i = float(lines[2][8:16])
raan = float(lines[2][17:25])
e = float("0." + lines[2][26:33])
argp = float(lines[2][34:42])
ma = float(lines[2][43:51])
n = float(lines[2][52:63])
t0 = lines[1][18:32]

# Print results to confirm what has been passed
print(f'i = {i}\nraan = {raan}\ne = {e}\nargp = {argp}\nma = {ma}\nn = {n}\n')

# Convert to appropriate units
i = deg_to_rad(i)      # Inclination (rad)
raan = deg_to_rad(raan)  # Longitude of Ascending Node (rad)
argp = deg_to_rad(argp)  # Argument of Perigee (rad)
ma = deg_to_rad(ma)   # Mean Anomaly (rad)
n = convert_n(n) # Mean Motion (rad/s)

a = compute_semimajoraxis(n,mu)
period = compute_period(mu,a)
print(f'period={period}')

# Create a time series (seconds) to model the XMM orbit on (uncomment for use)
# t = np.linspace(0,604800,10000)   # One week
# t = np.linspace(0,7200,10000)   # Two hours
# t = np.linspace(0,4000,10000)   # 80 minutes
n_period = 3
n_div = 1000
t = np.linspace(0,n_period*period,n_div)
dt = n_period*period/n_div

ta, a, peri, r, v_r, v_n, r_j2 = tle_orbit(t, mu, i,raan,e,argp,ma,n,dt)

# Calculate the altitude
alt = ot.magnitude(peri) - r_E

# Plot tangential velocity
plt.figure()
plt.title("XMM Tangential velocity")
plt.xlabel("Time (s)")
plt.ylabel("Tangential velocity (m/s)")
plt.plot(t,v_n,color='r')
plt.savefig("velocity.png")
plt.show()

# Plot altitude
plt.figure()
plt.title("XMM Altitude")
plt.xlabel("Time (s)")
plt.ylabel("Altitude (km)")
plt.plot(t,alt/1000.0,color='r')
plt.savefig("altitude.png")
plt.show()

# Plot on 3D graph

# Earth
theta = np.linspace(0, 2.*np.pi, 100)
phi = np.linspace(0, np.pi, 100)
x = r_E/1000000.0 * np.outer(np.cos(theta), np.sin(phi))
y = r_E/1000000.0 * np.outer(np.sin(theta), np.sin(phi))
z = r_E/1000000.0 * np.outer(np.ones(np.size(theta)), np.cos(phi))

r_ecef = ot.eci_to_ecef(t,r,t0)
r_j2_ecef = ot.eci_to_ecef(t,r_j2,t0)

# Perifocal Frame
plt.figure()
ax = plt.axes(projection='3d')
ax.set_aspect('equal')
set_aspect_equal(ax,peri[0]/1000000.0,peri[1]/1000000.0,peri[2]/1000000.0)
ax.plot_surface(x, y, z, color='g',alpha=0.1)
ax.scatter(peri[0]/1000000.0,peri[1]/1000000.0,peri[2]/1000000.0,color='b',marker='.',s=0.2)
ax.set_title("XMM Orbital Trace in Perifocal Frame")
ax.set_xlabel(r'$Position, p  (\times 10^3 km)$')
ax.set_ylabel(r'$Position, q  (\times 10^3 km)$')
ax.set_zlabel(r'$Position, w  (\times 10^3 km)$')
ax.ticklabel_format(style='plain')
plt.savefig("perifocal.png")
plt.show()

# ECI Frame
plt.figure()
ax = plt.axes(projection='3d')
ax.set_aspect('equal')
set_aspect_equal(ax,r[0]/1000000.0,r[1]/1000000.0,r[2]/1000000.0)
ax.plot_surface(x, y, z, color='g',alpha=0.1)
ax.scatter(r[0]/1000000.0,r[1]/1000000.0,r[2]/1000000.0,color='b',marker='.',s=0.2)
ax.set_title("XMM Orbital Trace in ECI Frame")
ax.set_xlabel(r'$Position, x  (\times 10^3 km)$')
ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
ax.ticklabel_format(style='plain')
plt.savefig("eci.png")
plt.show()

# ECEF Frame
plt.figure()
ax = plt.axes(projection='3d')
ax.set_aspect('equal')
set_aspect_equal(ax,r_ecef[0]/1000000.0,r_ecef[1]/1000000.0,r_ecef[2]/1000000.0)
ax.plot_surface(x, y, z, color='g',alpha=0.1)
ax.scatter(r_ecef[0]/1000000.0,r_ecef[1]/1000000.0,r_ecef[2]/1000000.0,c=t,marker='.',s=0.2)
ax.set_title("XMM Orbital Trace in ECEF Frame")
ax.set_xlabel(r'$Position, x  (\times 10^3 km)$')
ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
ax.ticklabel_format(style='plain')
plt.savefig("ecef.png")
plt.show()


# Convert to latitude and longitude
long = np.arctan2(r_ecef[1],r_ecef[0])/np.pi * 180
lat = np.arctan2(r_ecef[2],np.sqrt(r_ecef[0]**2+r_ecef[1]**2))/np.pi * 180

plt.figure()
plt.ylim(-90,90)
plt.xlim(-179,180)
plt.scatter(long,lat,color='b',marker='.',s=0.5)
plt.title("XMM Ground Track")
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.savefig("ground_track.png")
plt.show()

# J2 Perturbation
# ECEF Frame
plt.figure()
ax = plt.axes(projection='3d')
ax.set_aspect('equal')
set_aspect_equal(ax,r_j2_ecef[0]/1000000.0,r_j2_ecef[1]/1000000.0,r_j2_ecef[2]/1000000.0)
ax.plot_surface(x, y, z, color='g',alpha=0.1)
ax.scatter(r_j2_ecef[0]/1000000.0,r_j2_ecef[1]/1000000.0,r_j2_ecef[2]/1000000.0,c=t,marker='.',s=0.2)
ax.set_title("XMM Orbital Trace in ECEF Frame with J2 Perturbation")
ax.set_xlabel(r'$Position, x  (\times 10^3 km)$')
ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
ax.ticklabel_format(style='plain')
plt.savefig("j2_ecef.png")
plt.show()

# Convert to latitude and longitude
long_j2 = np.arctan2(r_j2_ecef[1],r_j2_ecef[0])/np.pi * 180
lat_j2 = np.arctan2(r_j2_ecef[2],np.sqrt(r_j2_ecef[0]**2+r_j2_ecef[1]**2))/np.pi * 180

plt.figure()
plt.ylim(-90,90)
plt.xlim(-179,180)
plt.scatter(long_j2,lat_j2,color='r',marker='.',s=0.1,label='With J2 Perturbation')
plt.scatter(long,lat,color='b',marker='.',s=0.1,label='No J2')
plt.title("XMM Ground Track with J2 Perturbation")
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.legend()
plt.savefig("ground_track_j2.png")
plt.show()

print("Q1 Complete")

########################################    QUESTION 2   ########################################

#### Execute Plane-change with Velocity-change manoeuvre as follows:
# 1. Elliptical transfer orbit on the equatorial plane from point on circular orbit to Ascending Node or Descending Node of the XMM orbit
# 2. Plane-change burn angle=inclination, including velocity-change (in one burn)

# Calculate position of Ascending Node (AN) and Descending Node (DN) of the XMM orbit
p_an, q_an, w_an, r_an,  p_dn, q_dn, w_dn, r_dn = compute_an_dn(i,raan,e,mu,a,argp)
# Distance of AN from Earth's centre
an_dist = ot.magnitude(r_an)
# Distance of DN from Earth's centre
dn_dist = ot.magnitude(r_dn)

print(f'AN at {an_dist[0]/1000.0}km and DN at {dn_dist[0]/1000.0}km')

### Try transfer orbit to the AN first

# Compute initial velocity
r_B = r_E + 500000
v_B1 = np.sqrt(mu/r_B)
# Compute angular momentum of transfer orbit to AN
r_C = an_dist
h_2 = np.sqrt(2*mu) * np.sqrt((r_B*r_C)/(r_B + r_C))
# Compute the speed required to enter transfer orbit
v_B2 = h_2/r_B
# Compute speed at apogee of transfer orbit/AN of target orbit
v_C2 = h_2/r_C

# Compute delta v to enter transfer orbit
dv_B = v_B2 - v_B1

# Compute speed at AN for the XMM orbit
v_r, v_n = ot.compute_orbital_velocity(a,e,-argp,mu)
v_C3 = v_n

# Compute delta v to plane-change from equatorial plane (0) to XMM plane (i) incl. velocity-change
dv_C = np.sqrt(v_C2**2 + v_C3**2 - 2*v_C2*v_C3*np.cos(i))

# Compute total delta v for plane-change at AN
dv_an = dv_B + dv_C
print(f'Delta v at AN = {dv_an/1000.0}km/s')

### Now try transfer orbit to the DN

# Compute initial velocity
r_B = r_E + 500000
v_B1 = np.sqrt(mu/r_B)
# Compute angular momentum of transfer orbit to DN
r_C = dn_dist
h_2 = np.sqrt(2*mu) * np.sqrt((r_B*r_C)/(r_B + r_C))
# Compute the speed required to enter transfer orbit
v_B2 = h_2/r_B
# Compute speed at apogee of transfer orbit/DN of target orbit
v_C2 = h_2/r_C

# Compute delta v to enter transfer orbit
dv_B = v_B2 - v_B1

# Compute speed at DN for the XMM orbit
v_r, v_n = ot.compute_orbital_velocity(a,e,np.pi-argp,mu)
v_C3 = v_n

# Compute delta v to plane-change from equatorial plane (0) to XMM plane (i) incl. velocity-change
dv_C = np.sqrt(v_C2**2 + v_C3**2 - 2*v_C2*v_C3*np.cos(i))

# Compute total delta v for plane-change at DN
dv_dn = dv_B + dv_C
print(f'Delta v at DN = {dv_B} + {dv_C} = {dv_dn/1000.0}km/s')

## Delta v at AN = [5.36414715]km/s
## Delta v at DN = [5.03422874]km/s
## Therefore we choose to plane-change at the DN for a more efficient transfer

# Points for the initial orbit to be plotted as a dotted line
angles = np.linspace(0,2*np.pi,10000)
x_init = r_B*np.cos(angles)
y_init = r_B*np.sin(angles)
z_init = np.zeros(len(angles))

# Compute transfer ellipse
ra2 = dn_dist   # Apogee
rp2 = r_B       # Perigee
a2 = (ra2+rp2)/2
e2 = compute_e(ra2,rp2)

# Compute specific angular momentum
h2 = compute_h(mu,ra2,rp2)

# Compute mean motion
n2 = compute_n(mu,a2)

# Determine the TLE to describe the transfer and final orbit
# Compute the orbits over
time = 3*24*60*60 # 3 days
divisions = 10000
dt = time/divisions
t = np.linspace(0,time,divisions)   
ta2, new_a2, peri2, r2, v_r2, v_n2, r_j2 = tle_orbit(t, mu, 0,raan-np.pi/2,e2,np.pi/2,0,n2) # starts from ma=0 (perigee), raan of transfer orbit is pi/2 from raan of final orbit, perigee is pi/2 from raan (aligned with raan of final orbit)
ta3, new_a3, peri3, r3, v_r3, v_n3, r_j2 = tle_orbit(t, mu, i,raan,e,argp,convert_ta_to_ma(argp,e),n) # starts from ta=argp (DN)

# Transpose to iterate by coordinate
peri2 = np.transpose(peri2)
peri3 = np.transpose(peri3)
r2 = np.transpose(r2)
r3 = np.transpose(r3)

# Take slices based on the true anomaly vectors

# Takes the transfer orbit up until the first apogee (ta ~= pi)
ta_current = -100
index = 0
peri_tran = np.array([peri2[0]])
r_tran = np.array([r2[0]])
v_tran_n = np.array([v_n2[0]])
v_tran_r = np.array([v_r2[0]])
while (ta_current > np.pi+0.01 or ta_current < np.pi-0.01):
    peri_tran = np.append(peri_tran,np.array([peri2[index]]),axis=0)
    r_tran = np.append(r_tran,np.array([r2[index]]),axis=0)
    v_tran_n = np.append(v_tran_n,v_n2[index])
    v_tran_r = np.append(v_tran_r,v_r2[index])
    ta_current = ta2[index]
    index += 1

### Takes the final orbit for 1 period
index = 0
peri_fin = np.array([peri3[0]])
r_fin = np.array([r3[0]])
v_fin_n = np.array(v_n3[0])
v_fin_r = np.array(v_r3[0])
while (t[index]<period):
    peri_fin = np.append(peri_fin,np.array([peri3[index]]),axis=0)
    r_fin = np.append(r_fin,np.array([r3[index]]),axis=0)
    v_fin_n = np.append(v_fin_n,v_n3[index])
    v_fin_r = np.append(v_fin_r,v_r3[index])
    index += 1

# Transpose back to 3 sets (p,q,w) or (x,y,z)
peri_tran = np.transpose(peri_tran)
r_tran = np.transpose(r_tran)
peri_fin = np.transpose(peri_fin)
r_fin = np.transpose(r_fin)

# Plot trajectory
plt.figure()
ax = plt.axes(projection='3d')
ax.set_aspect('equal')
set_aspect_equal(ax,r_fin[0]/1000000.0,r_fin[1]/1000000.0,r_fin[2]/1000000.0)
ax.plot_surface(x, y, z, color='g',alpha=0.1)
ax.scatter(x_init/1000000.0,y_init/1000000.0,z_init/1000000.0,color='g',marker='x',s=0.1,label='Initial orbit')
ax.scatter(r_tran[0]/1000000.0,r_tran[1]/1000000.0,r_tran[2]/1000000.0,color='b',marker='.',s=0.2,label='Transfer orbit')
ax.scatter(r_fin[0]/1000000.0,r_fin[1]/1000000.0,r_fin[2]/1000000.0,color='r',marker='.',s=0.2,label='Final orbit')
ax.legend(loc='upper right')
ax.set_title("XMM Orbital Transfer Trace in ECI Frame")
ax.set_xlabel(r'$Position, x  (\times 10^3 km)$')
ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
ax.ticklabel_format(style='plain')
plt.savefig("xmm_transfer_eci.png")
plt.show()

# Combine velocities into single vector
v_n = np.array([v_B1])
v_r = np.array([0])
t = np.array([0])

for i in range(1000):
    v_n = np.append(v_n,np.array([v_B1]))
    v_r = np.append(v_r,np.array([0]))
    t = np.append(t,np.array([t[-1]+dt]))

for i in range(len(v_tran_n)):
    v_n = np.append(v_n,np.array([v_tran_n[i]]))
    v_r = np.append(v_r,np.array([v_tran_r[i]]))
    t = np.append(t,np.array([t[-1]+dt]))

for i in range(len(v_fin_n)):
    v_n = np.append(v_n,np.array([v_fin_n[i]]))
    v_r = np.append(v_r,np.array([v_fin_r[i]]))
    t = np.append(t,np.array([t[-1]+dt]))

# Plot velocities over the entire transfer
plt.figure()
plt.title("Transfer to XMM Orbit velocity")
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.scatter(t,v_n,color='r',s=0.2,label='Tangential')
plt.scatter(t,v_r,color='b',s=0.2,label='Radial')
plt.legend(loc='upper right')
plt.savefig("transfer_velocity.png")
plt.show()

print('Q2 Complete')

########################################    QUESTION 3   ########################################

### Calculate and plot the Earth and Comet C/2023 E1's orbits around the Sun

## EARTH
i_earth = deg_to_rad(7.155)      # Inclination (rad)
raan_earth = deg_to_rad(174.9)  # Longitude of Ascending Node (rad)
argp_earth = deg_to_rad(288.1)  # Argument of Perigee (rad)
ma_earth = 0   # Mean Anomaly (rad)
e_earth = 0.0167086 # Eccentricity
a_earth = 149.60 * 10**9 # Semi-major axis (m)
n_earth = compute_n(mu_s,a_earth) # Mean Motion (rad/s)
period_earth = compute_period(mu_s,a_earth)

# Create a time series (seconds) to model the Earth orbit on
n_period = 1
n_div = 1000
t = np.linspace(0,n_period*period_earth,n_div)
dt = n_period*period/n_div

ta_earth, a_earth, peri_earth, r_earth, v_r_earth, v_n_earth, r_j2_earth = tle_orbit(t, mu_s, i_earth,raan_earth,e_earth,argp_earth,ma_earth,n_earth)

## COMET C/2023 E1
i_comet = deg_to_rad(38.31445216)      # Inclination (rad)
raan_comet = deg_to_rad(164.57526493741)  # Longitude of Ascending Node (rad)
argp_comet = deg_to_rad(105.89370862609)  # Argument of Perigee (rad)
ma_comet = 0   # Mean Anomaly (rad)
e_comet = 0.94693891 # Eccentricity
a_comet = 2894.366612 * 10**9 # Semi-major axis (m)
n_comet = compute_n(mu_s,a_comet) # Mean Motion (rad/s)
period_comet = compute_period(mu_s,a_comet)

# Create a time series (seconds) to model the comet orbit on
n_period = 1
n_div = 10000
t = np.linspace(0,n_period*period_comet,n_div)
dt = n_period*period/n_div

ta_comet, a_comet, peri_comet, r_comet, v_r_comet, v_n_comet, r_j2_comet = tle_orbit(t, mu_s, i_comet,raan_comet,e_comet,argp_comet,ma_comet,n_comet)

## Sun
theta = np.linspace(0, 2.*np.pi, 100)
phi = np.linspace(0, np.pi, 100)
x = r_S * np.outer(np.cos(theta), np.sin(phi))
y = r_S * np.outer(np.sin(theta), np.sin(phi))
z = r_S * np.outer(np.ones(np.size(theta)), np.cos(phi))

### Calculate Hohmann Transfer ellipse

# Calculate position of Ascending Node (AN) and Descending Node (DN) of COMET C/2023 E1 orbit
p_an, q_an, w_an, r_an,  p_dn, q_dn, w_dn, r_dn = compute_an_dn(i_comet,raan_comet,e_comet,mu_s,a_comet,argp_comet)
# Distance of AN from Earth's centre
an_dist = ot.magnitude(r_an)
# Distance of DN from Earth's centre
dn_dist = ot.magnitude(r_dn)
print(f'AN at {an_dist[0]/1000.0}km and DN at {dn_dist[0]/1000.0}km')

# AN is chosen as it is further away => more efficient transfer

# Compute angular momentum of transfer orbit to AN
r_B = a_earth # Can approximate Earth orbit as a circle with r_earth ~ a_earth
v_B1 = np.sqrt(mu_s/r_B)
r_C = an_dist[0]
h_2 = np.sqrt(2*mu_s) * np.sqrt((r_B*r_C)/(r_B + r_C))
# Compute the speed required to enter transfer orbit (equivalent to hyperbolic escape excess velocity)
v_excess = h_2/r_B - v_B1
# Compute speed at apogee of transfer orbit/AN of target orbit
v_C2 = h_2/r_C

# Compute speed at AN for the Comet C/2023 E1 orbit
v_r, v_n = ot.compute_orbital_velocity(a_comet,e_comet,-argp_comet,mu_s)
v_C3 = v_n

# Compute delta v to plane-change from equatorial plane (0) to XMM plane (i) incl. velocity-change
dv_2 = np.sqrt(v_C2**2 + v_C3**2 - 2*v_C2*v_C3*np.cos(i_comet))

print(f'Plane-change + velocity-change burn = {dv_2/1000.0}km/s')
print(f'Excess velocity required to enter Hohmann Transfer = {v_excess/1000.0}km/s')


### Given this excess velocity needed, calculate the hyperbolic escape from Earth's sphere of influence
r_p = 600000 + r_E # 600km altitude as specified
e_hyp = 1 + r_p*v_excess**2/mu # eccentricity of hyperbola
h_hyp = r_p * np.sqrt(v_excess**2 + (2*mu)/r_p) # angular momentum of hyperbolic trajectory
v_p = h_hyp/r_p # tangential velocity at periapsis of hyperbolic trajectory

v_c = np.sqrt(mu/r_p) # tangential velocity at circular parking orbit

dv_1 = v_p - v_c

dv = dv_1 + dv_2

# Compute hyperbola parameters
a_hyp = (h_hyp**2/mu) * (1/(e_hyp**2-1))
b_hyp = a_hyp * np.sqrt(e_hyp**2-1)
beta_hyp = np.arccos(1/e_hyp)

# Hyperbolic escape trajectory
theta_hyp = np.linspace(0,(11/17)*np.pi,1000)
x_hyp = -a_hyp * (e_hyp+np.cos(theta_hyp))/(1+e_hyp*np.cos(theta_hyp)) + r_p + a_hyp
y_hyp = b_hyp * (np.sqrt(e_hyp**2-1)*np.sin(theta_hyp))/(1+e_hyp*np.cos(theta_hyp))
r_hyp = np.row_stack((x_hyp,y_hyp))
r_hyp = np.transpose(r_hyp)

# Rotate such that escape direction parallel to Earth's heliocentric velocity vector (x-axis)
rot_matrix = np.array([[np.cos(np.pi+beta_hyp),-np.sin(np.pi+beta_hyp)],[np.sin(np.pi+beta_hyp),np.cos(np.pi+beta_hyp)]])
# rot_matrix = np.array([[np.cos(np.pi),-np.sin(np.pi)],[np.sin(np.pi),np.cos(np.pi)]])
for i in range(len(r_hyp)):
    r_hyp[i] = np.dot(rot_matrix,r_hyp[i])
r_hyp = np.transpose(r_hyp)

# Circular parking orbit
theta_park = np.linspace(0,2*np.pi,1000)
x_park = r_p * np.cos(theta_park)
y_park = r_p * np.sin(theta_park)

# Plot initial orbit and escape trajectory
# We define the positive x-axis of a non-ECI Earth-centred frame to be in the direction of Earth's heliocentric velocity vector
plt.figure()
plt.gca().set_aspect('equal')
plt.title("Initial orbit and Hyperbolic escape trajectory")
plt.xlabel("x (10^3 km)")
plt.ylabel("y (10^3 km)")
plt.scatter(r_hyp[0]/1000000.0,r_hyp[1]/1000000.0,color='r',s=0.2)
plt.scatter(x_park/1000000.0,y_park/1000000.0,color='b',s=0.2)
plt.savefig("escape.png")
plt.show()


# Compute transfer ellipse
ra2 = r_C   # Apoapsis
rp2 = r_B       # Periapsis
a2 = (ra2+rp2)/2
e2 = compute_e(ra2,rp2)
# Compute specific angular momentum
h2 = compute_h(mu_s,ra2,rp2)
# Compute mean motion
n2 = compute_n(mu_s,a2)
period2 = compute_period(mu_s,a2)
# Determine the TLE to describe the transfer and final orbit
# Compute the orbits over
divisions = 1000

time_tran = period2
dt_tran = time_tran/divisions
t_tran = np.linspace(0,time_tran,divisions)

time_comet = period_comet
dt_comet = time_comet/divisions
t_comet = np.linspace(0,time_comet,divisions)
ta2, new_a2, peri2, r2, v_r2, v_n2, r_j2 = tle_orbit(t_tran, mu_s, 0,raan_comet+np.pi/2,e2,np.pi/2,0,n2) # starts from ma=0 (perigee), raan of transfer orbit is pi/2 from raan of final orbit, perigee is pi/2 from raan (aligned with raan of final orbit)
ta3, new_a3, peri3, r3, v_r3, v_n3, r_j2 = tle_orbit(t_comet, mu_s, i_comet,raan_comet,e,argp,convert_ta_to_ma(-argp,e_comet),n_comet) # starts from ta=-argp (AN)
# Transpose to iterate by coordinate
peri2 = np.transpose(peri2)
r2 = np.transpose(r2)
peri3 = np.transpose(peri3)
r3 = np.transpose(r3)


# Take slices based on the true anomaly vectors
# Takes the transfer orbit up until the first apogee (ta ~= pi)
ta_current = -100
index = 0
peri_tran = np.array([peri2[0]])
r_tran = np.array([r2[0]])
v_tran_n = np.array([v_n2[0]])
v_tran_r = np.array([v_r2[0]])
while (ta_current > np.pi+0.01 or ta_current < np.pi-0.01):
    peri_tran = np.append(peri_tran,np.array([peri2[index]]),axis=0)
    r_tran = np.append(r_tran,np.array([r2[index]]),axis=0)
    v_tran_n = np.append(v_tran_n,v_n2[index])
    v_tran_r = np.append(v_tran_r,v_r2[index])
    ta_current = ta2[index]
    index += 1

# # Takes the final orbit for 1 period
index = 0
peri_fin = np.array([peri3[0]])
r_fin = np.array([r3[0]])
v_fin_n = np.array(v_n3[0])
v_fin_r = np.array(v_r3[0])
while (t_comet[index]<period_comet):
    peri_fin = np.append(peri_fin,np.array([peri3[index]]),axis=0)
    r_fin = np.append(r_fin,np.array([r3[index]]),axis=0)
    v_fin_n = np.append(v_fin_n,v_n3[index])
    v_fin_r = np.append(v_fin_r,v_r3[index])
    index += 1

# Transpose back to 3 sets (p,q,w) or (x,y,z)
peri_tran = np.transpose(peri_tran)
r_tran = np.transpose(r_tran)
peri_fin = np.transpose(peri_fin)
r_fin = np.transpose(r_fin)

# Plot Hohmann Transfer on 3D graph
# SCCEF Frame
plt.figure()
ax = plt.axes(projection='3d')
ax.set_aspect('equal')
set_aspect_equal(ax,r_comet[0]/1000000000.0/8,r_comet[1]/1000000000.0/8,r_comet[2]/1000000000.0/8) # Fit in an eighth of the comet orbit
ax.plot_surface(x/1000000000.0, y/1000000000.0, z/1000000000.0, color='orange',alpha=0.9)
ax.scatter(r_earth[0]/1000000000.0,r_earth[1]/1000000000.0,r_earth[2]/1000000000.0,color='g',marker='.',s=0.2)
ax.scatter(r_tran[0]/1000000000.0,r_tran[1]/1000000000.0,r_tran[2]/1000000000.0,color='b',marker='.',s=0.2)
ax.scatter(r_comet[0]/1000000000.0,r_comet[1]/1000000000.0,r_comet[2]/1000000000.0,color='r',marker='.',s=0.2)
ax.set_title("Orbital Trace in SCCEF Frame")
ax.set_xlabel(r'$Position, x  (\times 10^6 km)$')
ax.set_ylabel(r'$Position, y  (\times 10^6 km)$')
ax.set_zlabel(r'$Position, z  (\times 10^6 km)$')
ax.ticklabel_format(style='plain')
plt.savefig("sccef.png")
plt.show()

# Combine velocities into single vector
v_n = np.array([v_B1])
v_r = np.array([0])
t = np.array([0])

for i in range(len(v_tran_n)):
    v_n = np.append(v_n,np.array([v_tran_n[i]]))
    v_r = np.append(v_r,np.array([v_tran_r[i]]))
    t = np.append(t,np.array([t[-1]+dt_tran]))

for i in range(len(v_fin_n)):
    v_n = np.append(v_n,np.array([v_fin_n[i]]))
    v_r = np.append(v_r,np.array([v_fin_r[i]]))
    t = np.append(t,np.array([t[-1]+dt_comet]))

# Plot velocities over the entire transfer
plt.figure()
plt.title("Transfer to Comet C/2023 E1 Orbit velocity")
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.scatter(t,v_n,color='r',s=0.2,label='Tangential')
plt.scatter(t,v_r,color='b',s=0.2,label='Radial')
plt.legend(loc='upper right')
plt.savefig("comet_transfer_velocity.png")
plt.show()

print('Q3 Complete')