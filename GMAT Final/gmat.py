import FindStateVector as fsv
import numpy as np

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
    
np.set_printoptions(suppress=True)

mu_s = 1.327 * 10**20               # Standard gravitational parameter (Sun) (m^3/s^2)


print('***************** GMAT MISSION PART 1 *****************')
jd1 = fsv.date_to_JD(0,1,2,2041)
r_earth1, v_earth1, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd1,'Earth')
print(f'r_earth = {np.array(r_earth1)}, v_earth = {np.array(v_earth1)} at {jd1}')
#   -98478977.8  109677506 -10211.5276         -22.6483182  -20.01337474   0.00186335
# print(f'SMA = {a} km\nECC = {e}\nINC = {i/np.pi * 180.0} deg\nRAAN = {raan/np.pi * 180.0} deg\nAOP = {argp/np.pi * 180.0} deg\nTA = {ta/np.pi * 180.0} deg')

r_mars1, v_mars1, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd1,'Mars')
print(f'r_mars = {np.array(r_mars1)}, v_mars = {np.array(v_mars1)} at {jd1}')
#    80194745.4 -195857501 -6069478.13        23.33248642 11.26292024 -0.33533112
# print(f'SMA = {a} km\nECC = {e}\nINC = {i/np.pi * 180.0} deg\nRAAN = {raan/np.pi * 180.0} deg\nAOP = {argp/np.pi * 180.0} deg\nTA = {ta/np.pi * 180.0} deg')

r_starman1, v_starman1, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd1,'Starman')
print(f'r_starman = {np.array(r_starman1)}, v_starman = {np.array(v_starman1)} at {jd1}')
#   -157103526  3272049.94 -1954156.39          -5.34757516 -31.50207474  -0.50078033

jd2 = fsv.date_to_JD(0,4,5,2042)
r_earth2, v_earth2, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd2,'Earth')
r_mars2, v_mars2, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd2,'Mars')
r_starman2, v_starman2, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd2,'Starman')

v_dep_earth, v_arr_starman = lamberts_problem(r_earth1, r_starman2, (jd2 - jd1)*24*60*60 , mu_s/(10**9), "pro", 0.1)

print(f'v_dep_earth = {v_dep_earth}')

print(v_dep_earth-v_earth1)
print(np.linalg.norm(v_dep_earth-v_earth1))   

print(f'Starman insertion dv = {v_starman2-v_arr_starman}')
print(np.linalg.norm(v_starman2-v_arr_starman))   

# Subtract GMAT UTCModJulian offset
print(f'Arrive Starman: {jd2-2430000.0}')



jd3 = fsv.date_to_JD(0,8,5,2042)
print(f'jd3 = {jd3-2430000.0}')

jd4 = fsv.date_to_JD(0,17,3,2043)
print(f'jd4 = {jd4-2430000.0}')

r_starman3, v_starman3, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd3,'Starman')
r_mars4, v_mars4, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd4,'Mars')
print(f'r_mars4 = {r_mars4}')

v_dep_starman, v_arr_mars = lamberts_problem(r_starman3, r_mars4, (jd4 - jd3)*24*60*60 , mu_s/(10**9), "pro", 1)

print(f'Starman-Mars transfer insertion dv = {v_dep_starman-v_starman3}')
print(np.linalg.norm(v_dep_starman-v_starman3))

print(f'Mars excess v = {v_arr_mars-v_mars4}')
print(np.linalg.norm(v_arr_mars-v_mars4))   


v_dep_starman, v_arr_mars = lamberts_problem(r_starman2, r_mars4, (jd4 - jd2)*24*60*60 , mu_s/(10**9), "pro", -0.1)

print(f'Starman-Mars transfer insertion dv = {v_dep_starman-v_starman2}')
print(np.linalg.norm(v_dep_starman-v_starman2))

print(f'Mars excess v = {v_mars4-v_arr_mars}')
print(np.linalg.norm(v_mars4-v_arr_mars))   

print('***************** GMAT MISSION PART 2 *****************')
# Start at mars
jd1 = fsv.date_to_JD(0,15,7,2043)
r_earth1, v_earth1, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd1,'Earth')
print(f'r_earth = {np.array(r_earth1)}, v_earth = {np.array(v_earth1)} at {jd1}')
#   56587727.1 -141144492  13921.8955        27.16377729 10.97217722 -0.00108225
# print(f'SMA = {a} km\nECC = {e}\nINC = {i/np.pi * 180.0} deg\nRAAN = {raan/np.pi * 180.0} deg\nAOP = {argp/np.pi * 180.0} deg\nTA = {ta/np.pi * 180.0} deg')

r_mars1, v_mars1, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd1,'Mars')
print(f'r_mars = {np.array(r_mars1)}, v_mars = {np.array(v_mars1)} at {jd1}')
#    117177373 189835473 1110310.99        -19.70074614  14.79105334   0.79242948
# print(f'SMA = {a} km\nECC = {e}\nINC = {i/np.pi * 180.0} deg\nRAAN = {raan/np.pi * 180.0} deg\nAOP = {argp/np.pi * 180.0} deg\nTA = {ta/np.pi * 180.0} deg')

jd2 = fsv.date_to_JD(0,29,6,2044)
r_earth2, v_earth2, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd2,'Earth')
r_mars2, v_mars2, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd2,'Mars')

v_dep_mars, v_arr_earth = lamberts_problem(r_mars1, r_earth2, (jd2 - jd1)*24*60*60 , mu_s/(10**9), "pro", 0.1)

print(f'v_dep_mars = {v_dep_mars}')
# -17.0950608   14.54104162  -0.12305106

print(v_dep_mars-v_mars1)
print(np.linalg.norm(v_dep_mars-v_mars1))   

# Subtract GMAT UTCModJulian offset
print(f'Arrive Earth: {jd2-2430000.0}')