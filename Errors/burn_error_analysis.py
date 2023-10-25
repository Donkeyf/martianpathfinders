import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv
import PorkchopPlot as pp
import findOrbitalElements as foe
import orbitalTransforms2 as ot
import compute_params as cp

########################################    CONSTANTS DEFINITIONS   ########################################
mu_s = 1.327 * 10**20               # Standard gravitational parameter (Sun) (m^3/s^2)

def julian_day_number_to_gregorian(jdn):
        """Convert the Julian Day Number to the proleptic Gregorian Year, Month, Day."""
        L = jdn + 68569
        N = int(4 * L / 146_097)
        L = L - int((146097 * N + 3) / 4)
        I = int(4000 * (L + 1) / 1_461_001)
        L = L - int(1461 * I / 4) + 31
        J = int(80 * L / 2447)
        day = L - int(2447 * J / 80) + 1
        L = int(J / 11)
        month = J + 2 - 12 * L
        year = 100 * (N - 49) + I + L
        thirtydays = [9,4,6,11]
        thirtyonedays = [1,3,5,7,8,10,12]
        
        if month == 2:
            if year%4 == 0:
                if day > 29:
                    month = 3
                    day = day - 29
            else:
                if day > 28:
                    month = 3
                    day = day - 28
        elif month in thirtydays:
            if day > 30:
                month = month + 1
                day = day - 30
        elif month in thirtyonedays:
            if day > 31:
                month = month + 1
                day = day - 31

        if month == 13:
            month = 1
            year = year + 1

        return dt.date(year,month,day)
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
## Porkchop Plot for Earth to Starman
p_A = 'Earth'
p_B = 'Starman'
earth_to_starman = pp.PorkchopPlot(p_A,   p_B,  0,21,3,2032, 24,24,3,2032,      10,                 400,500,                    10,               "pro",         120)
earth_to_starman.get_min_dv()

r1, v1 = fsv.find_orbital_elements(earth_to_starman.min_jd_1,p_A)
r2, v2 = fsv.find_orbital_elements(earth_to_starman.min_jd_2,p_B)

percent_errors = np.arange(-0.01,0.01,0.0001)
v_errors = np.array([v1])


for p in percent_errors:
    v_errors = np.append(v_errors,[v1-p*v1/(np.linalg.norm(v1))],axis=0)


h1, e1, a1, T1, n1, i1, raan1, argp1, ta1, ma1 = foe.find_orbital_elements(r1,v1,mu_s)
h2, e2, a2, T2, n2, i2, raan2, argp2, ta2, ma2 = foe.find_orbital_elements(r2,v2,mu_s)

print(f'h1={h1},e1={e1},a1={a1},T1={T1}.n1={n1},i1={i1},raan1={raan1},argp1={argp1},ta1={ta1},ma1={ma1}')
print(f'h2={h2},e2={e2},a2={a2},T2={T2}.n2={n2},i2={i2},raan2={raan2},argp2={argp2},ta2={ta2},ma2={ma2}')
tof = (ma2-ma1)/n1 # time of flight
print(f'Time of transfer = {tof}')

# ta, a, peri, r, v_r, v_n, v = cp.tle_orbit(np.array([0,tof]),mu_s,i1,raan1,e1,argp1,ma1,n1)

# distances = np.array([])
# for v_error in v_errors:
#     h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r1,v_error,mu_s)
#     print(f'{mu_s},{i},{raan},{e},{argp},{ma},{n}')
#     ta, a, peri, r, v_r, v_n, v = cp.tle_orbit(np.array([tof]),mu_s,i,raan,e,argp,ma,n)
#     print(v)
#     r = np.transpose(r)
#     distances = np.append(distances,np.linalg.norm(r[0]-r2))
#     print(f'{r[0]}-{r2}')

# distances = np.array([])
# velocities = np.array([])
# for v_error in v_errors:
#     h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r1,v_error,mu_s)
#     ma2 = ma + n*tof
#     print(f'{ma}-{ma2}')
#     p, q, w, dp, dq, dw = ot.elements_to_perifocal([ta],a,e,mu_s)
#     peri = np.row_stack((p,q,w))
#     v_peri = np.row_stack((dp,dq,dw))

#     # Convert Perifocal Frame to ECI Frame
#     transform = ot.perifocal_to_eci_matrix(i,raan,argp)    
#     # Apply transform to create ECI frame coordinates array
#     r = np.dot(transform,peri)
#     r = np.transpose(r)
#     v = np.dot(transform,v_peri)
#     v = np.transpose(v)

#     distances = np.append(distances,np.linalg.norm(r[0]-r2))
#     velocities = np.append(velocities,np.linalg.norm(v[0]-v2))

distances = np.array([])
velocities = np.array([])
print(r1,v_errors[0])
for i in range(len(v_errors)):
    v_error = v_errors[i]
    h, e, a, T, n, i, raan, argp, ta, ma = foe.find_orbital_elements(r1,v_error,mu_s)
    if ma == None:
        np.delete(percent_errors,i)
        continue
    h = np.sqrt(a*mu_s*(1-e**2))

    ma2 = ma + n*tof
    # print(f'{ma}-{ma2}')
    p, q, w, dp, dq, dw = ot.elements_to_perifocal([ta],a,e,mu_s)
    peri = np.row_stack((p,q,w))
    v_peri = np.row_stack((dp,dq,dw))

    # Convert Perifocal Frame to ECI Frame
    transform = ot.perifocal_to_eci_matrix(i,raan,argp)    
    # Apply transform to create ECI frame coordinates array
    r = np.dot(transform,peri)
    r = np.transpose(r)
    v = np.dot(transform,v_peri)
    v = np.transpose(v)

    distances = np.append(distances,np.linalg.norm(r[0]-r2))
    velocities = np.append(velocities,np.linalg.norm(v[0]-v2))


plt.figure()
plt.title('Change in arrival distance with burn errors')
plt.ylabel('Distance error (km)')
plt.xlabel('Percentage error in delta-v (%)')
plt.plot(percent_errors,distances[0:-1])
plt.savefig('earth_to_starman_burn_error.png')
plt.show()

