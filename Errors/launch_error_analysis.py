import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv
import PorkchopPlot as pp

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

def timingError(p_A,p_B, ut1, d1, m1, y1, ut2, d2, m2, y2):
    t = np.arange(dt.datetime(y1,m1,d1,ut1) - dt.timedelta(hours=48), dt.datetime(y1,m1,d1,ut1) +  dt.timedelta(hours=48), dt.timedelta(hours=1)).astype(dt.datetime)
    t_diffs = np.arange(0-48,0+48,1)

    jd_1s = np.array([])
    dvs = np.array([])
    jd_2 = fsv.date_to_JD(ut2,d2,m2,y2)
    r_B_2, v_B_2, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd_2 ,p_B)

    for depart_date in t:
        jd_1 = fsv.date_to_JD(depart_date.hour,depart_date.day,depart_date.month,depart_date.year)
        jd_1s = np.append(jd_1s,jd_1)
        r_A_1, v_A_1, h, e, a, T, n, i, raan, argp, ta, M = fsv.find_state_vectors(jd_1,p_A)
        v_dep_A, v_arr_B = lamberts_problem(r_A_1, r_B_2, (jd_2 - jd_1)*24*60*60 , mu_s/(10**9), "pro", 1)
        dvs = np.append(dvs,np.linalg.norm(v_dep_A-v_A_1) + np.linalg.norm(v_arr_B-v_B_2))


    plt.figure()
    plt.title(f'Change in Delta-v within launch window for {p_A}-{p_B}')
    plt.ylabel('Delta-v (km/s)')
    plt.xlabel('Error in manoeuvre time (hours)')
    plt.plot(t_diffs,dvs)
    plt.savefig(f'{p_A}_to_{p_B}_launch_error.png')
    plt.show()

timingError('Earth','Starman',0,1,2,2041,0,4,5,2042)

timingError('Starman','Mars',0,8,5,2042,0,17,3,2043)

timingError('Mars','Earth',0,15,7,2043,0,29,6,2044)