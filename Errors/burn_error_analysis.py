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

########################################    CONSTANTS DEFINITIONS   ########################################
mu_s = 1.327 * 10**20               # Standard gravitational parameter (Sun) (m^3/s^2)
MU_SUN = 1.327 * 10**11 # (km^3/s^2)

## Porkchop Plot for Earth to Starman
p_A = 'Earth'
p_B = 'Starman'
earth_to_starman = pp.PorkchopPlot(p_A,   p_B,  0,21,3,2032, 24,24,3,2032,      10,                 400,500,                    10,               "pro",         120)
earth_to_starman.get_min_dv()

print(f'Min dv = {earth_to_starman.min_dv} km/s')
tof = (earth_to_starman.min_jd_2-earth_to_starman.min_jd_1)*(24*60*60) # time of flight
print(f'Time of transfer = {tof/(24*60*60)} days')

# time propagation vector
t = np.linspace(0,tof,1000)

r1, v1 = fsv.find_state_vectors(earth_to_starman.min_jd_1,p_A)
r2, v2 = fsv.find_state_vectors(earth_to_starman.min_jd_2,p_B)

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


r_earth, v_earth = fsv.find_state_vectors(earth_to_starman.min_jd_1,p_A)
h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r_earth,v_earth,MU_SUN)
ta, a, peri, r_earth, v_r, v_n, v = cp.tle_orbit(t,mu_s,i,raan,e,argp,ma,n)

r_starman, v_starman = fsv.find_state_vectors(earth_to_starman.min_jd_1,p_B)
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

