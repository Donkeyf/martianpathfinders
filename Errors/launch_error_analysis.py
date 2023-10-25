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
## Porkchop Plot for Earth to Starman
p_A = 'Earth'
p_B = 'Starman'
earth_to_starman = pp.PorkchopPlot(p_A,   p_B,  0,21,3,2032, 24,24,3,2032,      100,                 400,500,                    103,               "pro",         120)
earth_to_starman.get_plot()
earth_to_starman.get_min_dv()


t = np.arange(julian_day_number_to_gregorian(int(earth_to_starman.min_jd_1)) - dt.timedelta(hours=48), julian_day_number_to_gregorian(int(earth_to_starman.min_jd_1)) +  dt.timedelta(hours=48), dt.timedelta(hours=1)).astype(dt.datetime)

jd_1s = np.array([])
dvs = np.array([])
jd_2 = earth_to_starman.min_jd_2

for depart_date in t:
    jd_1 = fsv.date_to_JD(depart_date.hour,depart_date.day,depart_date.month,depart_date.year)
    jd_1s = np.append(jd_1s,jd_1)
    r_A_1, v_A_1 = fsv.find_orbital_elements(jd_1,p_A)
    r_B_2, v_B_2 = fsv.find_orbital_elements(jd_2 ,p_B)
    v_dep_A, v_arr_B = lamberts_problem(r_A_1, r_B_2, (jd_2 - jd_1)*24*60*60 , mu_s/(10**9), "pro", 1)
    dvs = np.append(dvs,np.linalg.norm(v_dep_A-v_A_1) + np.linalg.norm(v_arr_B-v_B_2))


plt.figure()
plt.title('Change in Delta-v within launch window')
plt.ylabel('Delta-v (km/s)')
plt.xlabel('Departure date (Julian days)')
plt.plot(jd_1s,dvs)
plt.savefig('earth_to_starman_launch_error.png')
plt.show()





for i in range(21,26):
    print(f'Delta-v with launch date {i}/3/2032 = {earth_to_starman.get_dv(2032,3,i,2033,6,29)}')


i_array = np.array([])

done_min = 0
done_max = 0

for i in range(len(earth_to_starman.depart_dates)):
    if earth_to_starman.depart_dates[i] == julian_day_number_to_gregorian(int(earth_to_starman.min_jd_1)):
        i_array = np.append(i_array,i)
    if earth_to_starman.jd_1s[i] == earth_to_starman.min_jd_1 and not done_min:
        i_min = i
        done_min = 1
    if earth_to_starman.jd_2s[i] == earth_to_starman.min_jd_2 and not done_max:
        i_max = i
        done_max = 1

print(f'Day {i_min} to Day {i_max}')

dvs = np.array([])
jd_1s = np.array([])
jd_diff = np.array([])

jd_1_last = 0
for i in range(len(earth_to_starman.jd_1s)):
    if earth_to_starman.jd_1s[i] != jd_1_last:
        jd_1_last = earth_to_starman.jd_1s[i]
        dvs = np.append(dvs,earth_to_starman.dvs[i+i_max-i_min])
        jd_1s = np.append(jd_1s,earth_to_starman.jd_1s[i+i_max-i_min])
        jd_diff = np.append(jd_diff,earth_to_starman.jd_2s[i+i_max-i_min]-earth_to_starman.jd_1s[i+i_max-i_min])

plt.figure()
plt.scatter(jd_1s,dvs)
plt.ylabel('Delta-v (km/s)')
plt.xlabel('Departure date (Julian days)')
plt.plot(jd_1s,jd_diff,c='red')
plt.savefig('earth_to_starman_launch_error2.png')
plt.show()


dvs = np.array([])
jd_1s = np.array([])
jd_2s = np.array([])

for i in i_array:
    dvs = np.append(dvs,earth_to_starman.dvs[int(i)])
    jd_1s = np.append(jd_1s,earth_to_starman.jd_1s[int(i)])
    jd_2s = np.append(jd_2s,earth_to_starman.jd_2s[int(i)])

print(np.shape(dvs))
print(np.shape(jd_1s))

plt.figure()
plt.scatter(jd_1s,dvs,c=jd_2s)
plt.ylabel('Delta-v (km/s)')
plt.xlabel('Departure date (Julian days)')
plt.savefig('earth_to_starman_launch_error_all.png')
plt.show()

# min_dvs = np.zeros(4)
# min_depart_dates = np.array([])
# n = -1
# min_dv = 1000
# for i in range(len(earth_to_starman.jd_1s)):
#     if earth_to_starman.depart_dates[i] in min_depart_dates:
#         if earth_to_starman.dvs[i] < min_dv:
#             min_dvs[n] = earth_to_starman.dvs[i]
#     else:
#         n = n + 1
#         min_depart_dates = np.append(min_depart_dates,earth_to_starman.depart_dates[i])
#         min_dv = 1000


# print(f'Minimum delta-v of {min_dvs}')