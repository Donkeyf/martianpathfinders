import numpy as np
import scipy as sp
import orbitalTransforms2 as ot
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import sys
import format_figures as ff
import compute_params as cp
import kepler_orbit as ko

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


j2_E = 1.08262668 * 10**(-3)         # Earth J2 Perturbation Constant

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

def earth_depart(v_tr,rot_esc_earth_deg):
    print('########################################    DEPARTURE FROM EARTH SOI   ########################################')
    # general parameters
    r_soi_E = cp.soi(R_E,m_E,m_S) # radius of Earth SOI (m)

    v_earth = np.sqrt(mu_s/R_E)
    # v_tr = v_earth + 1 # CHANGE THIS - transfer orbit velocity (relative to Sun)
    v_esc_earth = v_tr - v_earth # required excess escape velocity from Earth SOI (m/s)

    alt_park_earth = 1000000 # altitude of Earth parking orbit (m)
    radius_park_earth = r_E + alt_park_earth # radius of Earth parking orbit (m)

    print(f'r_soi_E = {r_soi_E/1000.0}km')
    
    earth_park = ko.KeplerOrbit('Earth',r_E,mu_e,radius_park_earth,22,0,0,-90,0,"31183.00000000",j2_E)
    n_period = 0.75
    n_div = 1000
    t = np.linspace(0,n_period*earth_park.period,n_div)
    dt = n_period*earth_park.period/n_div
    earth_park.tle_orbit(t,dt)

    exit(0)
    # SET parking orbit parameters HERE
    i_earth_park = cp.deg_to_rad(22)
    raan_earth_park = cp.deg_to_rad(0)
    e_earth_park = 0
    argp_earth_park = cp.deg_to_rad(-90)
    ma_earth_park = cp.deg_to_rad(0)
    n_earth_park = cp.compute_n(mu_e,radius_park_earth)
    t0_earth_park = "31183.00000000" # TODO: specific time of day

    a = cp.compute_semimajoraxis(n_earth_park,mu_e)
    period = cp.compute_period(mu_e,a)
    print(f'sma={a/1000.0}km')
    print(f'period={period}')

    # Create a time series (seconds) to model the XMM orbit on (uncomment for use)
    # t = np.linspace(0,604800,10000)   # One week
    # t = np.linspace(0,7200,10000)   # Two hours
    # t = np.linspace(0,4000,10000)   # 80 minutes
    n_period = 0.75
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
    print(f'e_hyp={e_hyp}\nbeta_hyp={beta_hyp}')

    # Hyperbolic escape trajectory in 3D
    theta_hyp = np.linspace(0,(14/17)*np.pi,1000)
    p_hyp = -a_hyp * (e_hyp+np.cos(theta_hyp))/(1+e_hyp*np.cos(theta_hyp)) + r_p + a_hyp
    q_hyp = b_hyp * (np.sqrt(e_hyp**2-1)*np.sin(theta_hyp))/(1+e_hyp*np.cos(theta_hyp))
    w_hyp = np.zeros(np.shape(q_hyp))


    peri_hyp = np.row_stack((p_hyp,q_hyp))
    peri_hyp = np.transpose(peri_hyp)
    # Rotate such that escape direction parallel to Earth's heliocentric velocity vector (x-axis)
    
    """
    For non-Hohmann transfer manoeuvre: hyperbolic escape direction is not parallel to earth velocity vector (ECI x-direction)
        - has inclination (CANNOT assume launch directly into parking orbit to escape  as above)
        - has rotation (implemented below)
    """
    rot_esc_earth_deg = 0
    rot_esc_earth = cp.deg_to_rad(rot_esc_earth_deg) # positive => inwards towards sun, negative => outwards away from sun
    rot_matrix = np.array([[np.cos(np.pi+rot_esc_earth-argp_earth_park),-np.sin(np.pi+rot_esc_earth-argp_earth_park)],[np.sin(np.pi+rot_esc_earth-argp_earth_park),np.cos(np.pi+beta_hyp+rot_esc_earth-argp_earth_park)]])
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
    cp.set_aspect_equal(ax,5*r[0]/1000000.0,5*r[1]/1000000.0,5*r[2]/1000000.0)
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

earth_depart(np.sqrt(mu_s/R_E)+1,0)