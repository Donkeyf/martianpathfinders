from math import nan
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv
import PorkchopPlot as pp
import findOrbitalElements as foe
import orbitalTransforms2 as ot
import compute_params as cp
import kepler_orbit as ko
import orbitalTransforms as ot
from delta_v_error import delta_v_error

########################################    CONSTANTS DEFINITIONS   ########################################
mu_e = 3.986004418 * 10**14          # Standard gravitational parameter (Earth) (m^3/s^2)
mu_s = 1.327 * 10**20               # Standard gravitational parameter (Sun) (m^3/s^2)
mu_m = 4.282837 * 10**13            # Standard gravitational parameter (Mars) (m^3/s^2)
MU_SUN = 1.327 * 10**11 # (km^3/s^2)


def lamberts_problem(R1, R2, t, mu, case="pro", z=0):
        # z: initial iteration value of solution for z


        # Stumpff functions
        def S(z):
            if z > 0:
                s = (np.sqrt(z) - np.sin(np.sqrt(z))) / np.sqrt(z)**3
            elif z < 0:
                s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / np.sqrt(-z)**3
            else:
                s = 1/6
            return s

        def C(z):
            if z > 0:
                c = (1 - np.cos(np.sqrt(z)))/z
            elif z < 0:
                c = (np.cosh(np.sqrt(-z)) - 1)/(-z)
            else:
                c = 1/2
            return c    

        # subfunctions
        def y(z):
            return r1 + r2 + A*(z*S(z) - 1)/np.sqrt(C(z))
        
        def F(z,t):
            F = (y(z)/C(z))**(1.5)*S(z) + A*np.sqrt(y(z)) - np.sqrt(mu)*t
            return F
        
        def dF(z):
            if z == 0:
                dF = np.sqrt(2)/40*y(0)**1.5 + A/8*( np.sqrt(y(0)) + A*np.sqrt(1/(2*y(0))) )
            else:
                dF = (y(z)/C(z))**1.5 * ( 1/2/z * (C(z)-3/2*S(z)/C(z)) + 3/4*S(z)**2/C(z) ) + A/8*( 3*S(z)/C(z)*np.sqrt(y(z)) + A*np.sqrt(C(z)/y(z)) )
            return dF


        # (1) calculate position magnitudes
        r1 = np.linalg.norm(R1)
        r2 = np.linalg.norm(R2)

        # (2) calculate change in true anomaly
        theta = np.arccos(np.dot(R1, R2) / r1 / r2)

        # z-component of positions cross product
        c_z = np.cross(R1, R2)[2]

        # check and update theta's quadrant as necessary
        if case == "pro":
            if c_z < 0:
                theta = 2*np.pi - theta
        elif case == "retro":
            if c_z >= 0:
                theta = 2*np.pi - theta

        # (3) 'A' constant
        A = np.sin(theta) * np.sqrt(r1*r2/(1-np.cos(theta)))

        # (4) solving for z with Newton's Method
        ## initialise z value
        while F(z,t) < 0:
            z += 0.1
            z = round(z, 5)

        ## Newton's method
        tol = 1                             # initialise tolerance value  
        while tol > 1e-8:                   # iterate until tolerance condition (<1e-8) is satisfied
            z_np1 = z - F(z,t)/dF(z)
            tol = abs(z_np1-z)              # update tolerance
            z = z_np1                       # update z for next iteration or for final answer

        # (5) Lagrange functions
        f = 1 - y(z)/r1
        g = A*np.sqrt(y(z)/mu)
        g_dot = 1 - y(z)/r2

        # (6) endpoint velocities of the transfer orbit
        V1 = 1/g * (R2 - f*R1)
        V2 = 1/g * (g_dot*R2 - R1)

        return V1, V2

def helio_to_eci(Vec_Helio, planet="Earth"):
    # convert vectors between Heliocentric frame to ECI frame
    # only applies for velocities, not positions 

    # obliquity (note negative angle for rotation)
    if planet == "Earth":
        o = -23.44 * np.pi/180
    elif planet == "Mars":
        o = -25.19 * np.pi/180

    # matrix rotation about the x-axis
    R_mat = np.matrix([[1, 0, 0], [0, np.cos(o), np.sin(o)], [0, -np.sin(o), np.cos(o)]])
    Vec_ECI = np.array(np.matmul(R_mat, Vec_Helio))[0]

    return Vec_ECI

def eci_to_helio(Vec_eci, planet="Earth"):
    # convert vectors between Heliocentric frame to ECI frame
    # only applies for velocities, not positions 

    # obliquity (note negative angle for rotation)
    if planet == "Earth":
        o = -23.44 * np.pi/180
    elif planet == "Mars":
        o = -25.19 * np.pi/180

    # matrix rotation about the x-axis
    R_mat_transpose = np.transpose(np.matrix([[1, 0, 0], [0, np.cos(o), np.sin(o)], [0, -np.sin(o), np.cos(o)]]))
    Vec_Helio = np.array(np.matmul(R_mat_transpose, Vec_eci))[0]

    return Vec_Helio

def burnError(p_A,p_B, mu, ut1, d1, m1, y1, ut2, d2, m2, y2, r_p_A, v_hyp_A, v_park_A):
    jd_1 = fsv.date_to_JD(ut1,d1,m1,y1) # departure date
    jd_2 = fsv.date_to_JD(ut2,d2,m2,y2) # arrival date

    r_A_1, v_A_1, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd_1,p_A) # departure planet at depart date
    r_B_2, v_B_2, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd_2,p_B) # arrival planet at arrive date
    
    # Lambert's solver
    tof = (jd_2 - jd_1)*(24*60*60) # time of flight
    print(f'Time of transfer = {tof/(24*60*60)} days')
    v_dep_A, v_arr_B = lamberts_problem(r_A_1, r_B_2, tof , mu_s/(10**9), "pro", 1)
    v_inf = v_dep_A - v_A_1
    print(f'v_inf = {v_inf} km/s')
    
    print(f'burn_dv = {v_hyp_A - v_park_A} km/s')
    v_inf2 = delta_v_error(r_p_A*1000.0,0,v_hyp_A*1000.0,mu)
    v_inf2 = eci_to_helio(v_inf2,p_A)
    print(f'alternative v_inf = {v_inf2/1000.0} km/s')

    # time propagation vector
    t = np.linspace(0,tof,1000)

    # Propagate two celestial bodies
    r_A, v_A, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd_1,p_A)
    h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r_A,v_A,MU_SUN)
    ta, a, peri, r_A_vec, v_r, v_n, v = cp.tle_orbit(t,mu_s,i,raan,e,argp,ma,n)
    r_B, v_B, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd_1,p_B)
    h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r_B,v_B,MU_SUN)
    ta, a, peri, r_B_vec, v_r, v_n, v = cp.tle_orbit(t,mu_s,i,raan,e,argp,ma,n)

    # Propagate target transfer
    h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r_A_1,v_dep_A,MU_SUN)
    ta, a, peri, r_transfer, v_r, v_n, v = cp.tle_orbit(t,mu_s,i,raan,e,argp,ma,n)

    # Final distance between target transfer and target planet
    r_transfer_last = np.array((r_transfer[0][-1],r_transfer[1][-1],r_transfer[2][-1]))
    r_B_vec_last = np.array((r_B_vec[0][-1],r_B_vec[1][-1],r_B_vec[2][-1]))
    print(f'r_transfer_last = {r_transfer_last}\nr_B_vec_last = {r_B_vec_last}\nr_B_2 = {r_B_2}')
    d = r_transfer_last - r_B_vec_last
    print(f'd = {d} km')

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_aspect('equal')
    cp.set_aspect_equal(ax,r_A_vec[0]/1000000000.0,r_A_vec[1]/1000000000.0,r_A_vec[2]/1000000000.0)
    ax.scatter(r_A_vec[0]/1000000000.0,r_A_vec[1]/1000000000.0,r_A_vec[2]/1000000000.0,color='b',marker='.',s=0.8)
    ax.scatter(r_B_vec[0]/1000000000.0,r_B_vec[1]/1000000000.0,r_B_vec[2]/1000000000.0,color='r',marker='.',s=0.8)
    ax.scatter(r_transfer[0]/1000000000.0,r_transfer[1]/1000000000.0,r_transfer[2]/1000000000.0,color='green',marker='.',s=0.8)
    ax.set_title(f"{p_A}-{p_B} Transfer Orbital Trace")
    ax.set_xlabel(r'$Position, x  (\times 10^6 km)$')
    ax.set_ylabel(r'$Position, y  (\times 10^6 km)$')
    ax.set_zlabel(r'$Position, z  (\times 10^6 km)$')
    ax.ticklabel_format(style='plain')
    plt.savefig(f'{p_A}_to_{p_B}_transfer.png')
    plt.show()

    percent_errors = np.linspace(-0.001,0.001,51) # -0.1% to +0.1% error
    d_vec = np.array([])
    for percent_error in percent_errors:
        v_inf = delta_v_error(r_p_A*1000.0,percent_error,v_hyp_A*1000.0,mu)
        print(f'v_inf = {v_inf}')
        v_dep_A = v_inf - v_A_1 
        print(f'r = {r_A_1}, v = {v_dep_A}')
        h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r_A_1/1000.0,v_dep_A/1000.0,MU_SUN)
        print(f'n, i, raan, argp, ta, ma = {n}, {i}, {raan}, {argp}, {ta}, {ma}')
        ta, a, peri, r_transfer, v_r, v_n, v = cp.tle_orbit(t,mu_s,i,raan,e,argp,ma,n)
        # Final distance between target transfer and target planet
        r_transfer_last = np.array((r_transfer[0][-1],r_transfer[1][-1],r_transfer[2][-1]))
        r_B_vec_last = np.array((r_B_vec[0][-1],r_B_vec[1][-1],r_B_vec[2][-1]))
        print(f'r_transfer_last = {r_transfer_last}\nr_B_vec_last = {r_B_vec_last}\nr_B_2 = {r_B_2}')
        d = r_transfer_last - r_B_vec_last
        d_vec = np.append(d_vec,np.linalg.norm(d))

    plt.figure()
    plt.title(f'Change in arrival distance with velocity errors for {p_A}-{p_B}')
    plt.ylabel('Distance (m)')
    plt.xlabel('Error in velocity after burn (%)')
    plt.plot(percent_errors*100.0,d_vec)
    plt.savefig(f'{p_A}_to_{p_B}_burn_error.png')
    plt.show()

burnError('Earth','Starman',mu_e,0,1,2,2041,0,4,5,2042,np.array([-429.0935112062964, 5835.709017133281, 3801.692126670696]),np.array([-11.061888968306045, 0.12664288989179331, -1.442945838248612]),np.array([-7.49396089161757, 0.08579518984232422, -0.977534642730493]))

burnError('Mars','Earth',mu_m,0,15,7,2043,0,29,6,2044,np.array([-1939.5777054982793, -2994.9653058026825, 1277.6232392111026]),np.array([4.3854639187425395, -3.2102280829675935, -0.8676843478721992]),np.array([2.6785819243801807, -1.9607638041268434, -0.5299698397574534]))

exit(0)

percent_errors = np.linspace(-0.001,0.001,5)
burn_dv = earth_to_starman.min_dv*v1/np.linalg.norm(v1)
v_errors = np.array([burn_dv])
for p in percent_errors:
    v_errors = np.append(v_errors,[v1+burn_dv+p*burn_dv],axis=0)
v_errors = np.delete(v_errors,0,axis=0)

v_errors_mag = np.array([])
for v_error in v_errors:
    v_errors_mag = np.append(v_errors_mag,np.linalg.norm(v_error))

# plt.figure()
# plt.title('Velocities')
# plt.scatter(percent_errors,v_errors_mag)
# plt.show()

d = np.zeros(len(v_errors))

r_vectors = np.zeros((len(v_errors),3,len(t)))

# Propagate each v_error over tof time series
for j in range(len(v_errors)):
    print(f'velcotiy = {v1+v_errors[j]}')
    h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r1,v_errors[j],MU_SUN)
    print(f'e = {e}')
    print('mu_s,i,raan,e,argp,ma,n')
    print(f'{mu_s},{i},{raan},{e},{argp},{ma},{n}')

    ta, a, peri, r_vectors[j], v_r, v_n, v = cp.tle_orbit(t,mu_s,i,raan,e,argp,ma,n)
    r_end = np.array([r_vectors[j][0][-1],r_vectors[j][1][-1],r_vectors[j][2][-1]])
    print(f'{r_end} - {r2}')
    d[j] = np.linalg.norm(r_end-r2)


r_earth, v_earth, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(earth_to_starman.min_jd_1,p_A)
h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r_earth,v_earth,MU_SUN)
ta, a, peri, r_earth, v_r, v_n, v = cp.tle_orbit(t,mu_s,i,raan,e,argp,ma,n)

r_starman, v_starman, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(earth_to_starman.min_jd_1,p_B)
h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r_starman,v_starman,MU_SUN)
ta, a, peri, r_starman, v_r, v_n, v = cp.tle_orbit(t,mu_s,i,raan,e,argp,ma,n)

# Plot
plt.figure()
ax = plt.axes(projection='3d')
ax.set_aspect('equal')
cp.set_aspect_equal(ax,r_vectors[0][0]/1000000.0,r_vectors[0][1]/1000000.0,r_vectors[0][2]/1000000.0)
ax.scatter(r_vectors[0][0]/1000000.0,r_vectors[0][1]/1000000.0,r_vectors[0][2]/1000000.0,color='r',marker='.',s=0.2)
ax.scatter(r_vectors[1][0]/1000000.0,r_vectors[1][1]/1000000.0,r_vectors[1][2]/1000000.0,color='b',marker='.',s=0.2)
ax.scatter(r_vectors[2][0]/1000000.0,r_vectors[2][1]/1000000.0,r_vectors[2][2]/1000000.0,color='green',marker='.',s=0.2)
ax.scatter(r_vectors[3][0]/1000000.0,r_vectors[3][1]/1000000.0,r_vectors[3][2]/1000000.0,color='orange',marker='.',s=0.2)
ax.scatter(r_vectors[4][0]/1000000.0,r_vectors[4][1]/1000000.0,r_vectors[4][2]/1000000.0,color='yellow',marker='.',s=0.2)
ax.scatter(r_earth[0]/1000000.0,r_earth[1]/1000000.0,r_earth[2]/1000000.0,color='b',marker='.',s=0.8)
ax.scatter(r_starman[0]/1000000.0,r_starman[1]/1000000.0,r_starman[2]/1000000.0,color='r',marker='.',s=0.8)
ax.set_title(f"Burn Error Orbital Trace in ECI Frame")
ax.set_xlabel(r'$Position, x  (\times 10^3 km)$')
ax.set_ylabel(r'$Position, y  (\times 10^3 km)$')
ax.set_zlabel(r'$Position, z  (\times 10^3 km)$')
ax.ticklabel_format(style='plain')
plt.savefig(f"burn_error_eci.png")
plt.show()

plt.figure()
plt.title('Change in arrival distance with burn errors')
plt.ylabel('Distance error (km)')
plt.xlabel('Percentage error in delta-v (%)')
plt.scatter(percent_errors*100,d)
plt.savefig('earth_to_starman_burn_error.png')
plt.show()

