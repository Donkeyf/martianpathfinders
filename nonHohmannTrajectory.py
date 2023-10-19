import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv


########################################    CONSTANTS DEFINITIONS   ########################################
mu_e = 3.986004418 * 10**14          # Standard gravitational parameter (Earth) (m^3/s^2)
mu_s = 1.327 * 10**20               # Standard gravitational parameter (Sun) (m^3/s^2)
mu_m = 4.282837 * 10**13            # Standard gravitational parameter (Mars) (m^3/s^2)

r_E = 6378000                      # Radius of Earth (m)
m_E = 5.974 * 10**24               # Mass of Earth (kg)
R_E = 149.6 * 10**9                 # Earth-Sun distance (m)

r_S = 696340000                    # Radius of the Sun (m)
m_S = 1.989 * 10**30               # Mass of the Sun (kg)

r_M = 3397000                      # Radius of Mars (m)
m_M = 5.974 * 10**24                # Mass of Mars (kg)
R_M = 227.9 * 10**9                 # Mars-Sun distance (m)


j2_E = 1.08262668 * 10**(-3)         # J2 Perturbation Constant for Earth
j2_M = 1.9555 * 10**(-3)               # J2 Perturbation Constant for Mars

########################################    Solve Lambert's problem     ########################################
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

# Testing the function with textbook example 5.2
mu_E = 3.986e5

R1 = np.array([5000, 10000, 2100])
R2 = np.array([-14600, 2500, 7000])

print('Example 5.2')
V1, V2 = lamberts_problem(R1, R2, 60*60, mu_E, "pro", 1)
print(V1, V2)

print()


# Depart Earth
jd_1 = fsv.date_to_JD(0,1,1,2031) # Change this date
r_earth_1, v_earth_1 = fsv.find_orbital_elements(jd_1,'Earth')
r_mars_1, v_mars_1 = fsv.find_orbital_elements(jd_1,'Mars')

# Arrive at Mars
jd_2 = fsv.date_to_JD(0,1,10,2032) # Change this date
r_earth_2, v_earth_2 = fsv.find_orbital_elements(jd_2,'Earth')
r_mars_2, v_mars_2 = fsv.find_orbital_elements(jd_2,'Mars')


print(f'Lambert\'s: {r_earth_1}, {r_mars_2}, {(jd_2 - jd_1)*24*60*60} , {mu_s/(10**9)}')

v_dep_earth, v_arr_mars = lamberts_problem(r_earth_1, r_mars_2, (jd_2 - jd_1)*24*60*60 , mu_s/(10**9), "pro", 1)



# Depart Starman
jd_3 = jd_2 + 4 * 30 * 24 * 60 * 60
# jd_3 = fsv.date_to_JD(0,7,2,33) # Change this date
r_earth_3, v_earth_3 = fsv.find_orbital_elements(jd_3,'Earth')
r_mars_3, v_mars_3 = fsv.find_orbital_elements(jd_3,'Mars')

# Arrive at Earth
jd_4 = fsv.date_to_JD(0,7,2,33) # Change this date
r_earth_4, v_earth_4 = fsv.find_orbital_elements(jd_4,'Earth')
r_mars_4, v_mars_4 = fsv.find_orbital_elements(jd_4,'Mars')

v_dep_mars, v_arr_earth = lamberts_problem(r_earth_3, r_mars_4, jd_2 - jd_1 , mu_s, "pro", 1)



print(f'Earth escape excess velocity = {v_dep_earth} km/s, Earth velocity = {v_earth_1}')
print(f'Mars capture excess velocity = {v_arr_mars} km/s, Mars velocity = {v_mars_2}')
print(f'Mars escape excess velocity = {v_dep_mars} km/s, Mars velocity = {v_mars_3}')
print(f'Earth capture excess velocity = {v_arr_earth} km/s, Earth velocity = {v_earth_4}')
print()
print(f'Earth escape excess velocity = {v_dep_earth-v_earth_1} km/s')
print(f'Mars capture excess velocity = {v_arr_mars-v_mars_2} km/s')
print(f'Mars escape excess velocity = {v_dep_mars-v_mars_3} km/s')
print(f'Earth capture excess velocity = {v_arr_earth-v_earth_4} km/s')
print()

print(np.linalg.norm(v_dep_earth-v_earth_1))
print(np.linalg.norm(v_arr_mars-v_mars_2))
print(np.linalg.norm(v_dep_mars-v_mars_3))
print(np.linalg.norm(v_arr_earth-v_earth_4))

plt.figure()
plt.scatter(r_earth_1[0],r_earth_1[1],color='b',alpha=0.2)
plt.scatter(r_earth_2[0],r_earth_2[1],color='b',alpha=0.4)
plt.scatter(r_earth_3[0],r_earth_3[1],color='b',alpha=0.6)
plt.scatter(r_earth_4[0],r_earth_4[1],color='b',alpha=0.8)

plt.scatter(r_mars_1[0],r_mars_1[1],color='orange',alpha=0.2)
plt.scatter(r_mars_2[0],r_mars_2[1],color='orange',alpha=0.4)
plt.scatter(r_mars_3[0],r_mars_3[1],color='orange',alpha=0.6)
plt.scatter(r_mars_4[0],r_mars_4[1],color='orange',alpha=0.8)
plt.show()


# Pork chop plot

jd_start = fsv.date_to_JD(0,1,1,2031) # Change this date
print(f'start date should be 1/1/2031: {julian_day_number_to_gregorian(int(jd_start))}')
jd_end = fsv.date_to_JD(0,1,10,2032) # Change this date
print(f'end date should be 1/10/2032: {julian_day_number_to_gregorian(int(jd_end))}')
print(f'{jd_start} - {jd_end}')
dates = np.linspace(jd_start,jd_end,60)
times = np.linspace(0.5*30,36*30,30)  # days

jd_1s = np.array([])
jd_2s = np.array([])
dvs = np.array([])

for depart in range(len(dates)):
    jd_1 = dates[depart]
    for time in times:
        jd_1s = np.append(jd_1s,jd_1)
        jd_2 = jd_1 + time
        jd_2s = np.append(jd_2s,jd_2)
        r_earth_1, v_earth_1 = fsv.find_orbital_elements(jd_1,'Earth')
        r_mars_2, v_mars_2 = fsv.find_orbital_elements(jd_2,'Mars')
        v_dep_earth, v_arr_mars = lamberts_problem(r_earth_1, r_mars_2, (jd_2 - jd_1)*24*60*60 , mu_s/(10**9), "pro", 1)
        dvs = np.append(dvs,np.linalg.norm(v_dep_earth-v_earth_1) + np.linalg.norm(v_arr_mars-v_mars_2))
        # print(f'dv = {np.linalg.norm(v_dep_earth-v_earth_1) + np.linalg.norm(v_arr_mars-v_mars_2)}')

print(f'dv range = {np.min(dvs)/1000.0} to {np.max(dvs)/1000.0} km/s')



# plt.figure()
# ax = plt.axes(projection='3d')
# ax.set_xlabel('Earth Departure Date')
# ax.set_ylabel('Mars Arrival Date')
# ax.set_zlabel('Delta-v (km/s)')
# ax.scatter(jd_1s/(60.0*60*24),jd_2s/(60.0*60*24),dvs,c=np.log(dvs))
# plt.show()

# Remove outliers
new_jd_1s = np.array([])
new_jd_2s = np.array([])
new_dvs = np.array([])

min_dv = 100

for i in range(len(jd_1s)):
    if dvs[i] < min_dv:
        min_dv = dvs[i]
        min_jd_1 = jd_1s[i]
        min_jd_2 = jd_2s[i]
    if dvs[i] < 120:
        new_jd_1s = np.append(new_jd_1s,jd_1s[i])
        new_jd_2s = np.append(new_jd_2s,jd_2s[i])
        new_dvs = np.append(new_dvs,dvs[i])

# Convert Julian times to DD/MM/YY

depart_dates = [julian_day_number_to_gregorian(int(new_jd_1s[i])) for i in range(len(new_jd_1s))]
arrive_dates = [julian_day_number_to_gregorian(int(new_jd_2s[i])) for i in range(len(new_jd_1s))]


# for i in range(len(new_jd_1s)):
#     depart_dates.append(dt.date(julian_day_number_to_gregorian(int(new_jd_1s[i]))))
#     arrive_dates.append(dt.date(julian_day_number_to_gregorian(int(new_jd_2s[i]))))


print(f'Minimum delta-v of {min_dv} km/s departing at {julian_day_number_to_gregorian(int(min_jd_1))} and arriving at {julian_day_number_to_gregorian(int(min_jd_2))}')


zeros = np.zeros(len(new_dvs))
plt.figure()
plt.title('Colour scale')
plt.ylabel('Delta-v (km/s)')
plt.scatter(zeros,new_dvs,c=np.log(new_dvs),s=50)
plt.savefig('colour_scale_eath_mars.png')
plt.show()

plt.figure()
plt.title('Porkchop plot for Earth-Mars Transfer')
plt.xlabel('Departure from Earth date')
plt.ylabel('Arrival to Mars date')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=40))
plt.gca().yaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
plt.gca().yaxis.set_major_locator(mdates.DayLocator(interval=60))
plt.scatter(depart_dates,arrive_dates,c=np.log(new_dvs),s=100)
plt.gcf().autofmt_xdate()
plt.savefig('porkchop_earth_mars.png')
plt.show()

# Get Delta-v of particular departure and arrival date
departure_get = dt.date(2031,7,1)
arrival_get = dt.date(2032,10,1)

closest_depart = dt.date(1000,1,1)
for depart in depart_dates:
    if abs(departure_get-depart) < abs(departure_get-closest_depart):
        closest_depart = depart

closest_arrive = dt.date(1000,1,1)
for i in range(len(arrive_dates)):
    if depart_dates[i] == closest_depart and abs(arrival_get-arrive_dates[i]) < abs(arrival_get-closest_arrive):
        closest_arrive = arrive_dates[i]

print(f'Closest dates: {closest_depart} - {closest_arrive}')

n = 0
while depart_dates[n] != closest_depart:
    n = n+1
while arrive_dates[n] != closest_arrive:
    n = n+1

print(f'Delta-v = {new_dvs[n]}')