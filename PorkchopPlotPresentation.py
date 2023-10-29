import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv


########################################    CONSTANTS DEFINITIONS   ########################################
mu_s = 1.327 * 10**20               # Standard gravitational parameter (Sun) (m^3/s^2)

class PorkchopPlot:
    def __init__(self, p_A, p_B, ut1, d1, m1, y1, ut2, d2, m2, y2, n_jds, t_start, t_end, n_ts, dir, dv_max):
        self.p_A = p_A
        self.p_B = p_B
        self.jd_start = fsv.date_to_JD(ut1,d1,m1,y1)
        self.jd_end = fsv.date_to_JD(ut2,d2,m2,y2)
        self.n_jds = n_jds
        self.t_start = t_start
        self.t_end = t_end
        self.n_ts = n_ts
        self.dir = dir
        self.dv_max = dv_max
        print(f'### Initializing Porkchop Plot for {self.p_A} to {self.p_B} departing at {PorkchopPlot.julian_day_number_to_gregorian(int(self.jd_start))} to {PorkchopPlot.julian_day_number_to_gregorian(int(self.jd_end))} with a flight time of {self.t_start}-{self.t_end} days, generating {self.n_jds*self.n_ts} orbits...')

        # Input time series: start dates and flight times
        dates = np.linspace(self.jd_start,self.jd_end,n_jds)
        times = np.linspace(t_start,t_end,n_ts)  # days

        # Output time series: departure and arrival dates
        self.jd_1s = np.array([])
        self.jd_2s = np.array([])
        # Output delta-v series
        self.dvs = np.array([])
        self.r_A = np.array([])
        self.r_B = np.array([])

        i = 0
        for depart in range(len(dates)):
            i += 1
            jd_1 = dates[depart]
            print(str(i) + "/" + str(len(dates)))
            for time in times:
                self.jd_1s = np.append(self.jd_1s,jd_1)
                jd_2 = jd_1 + time
                self.jd_2s = np.append(self.jd_2s,jd_2)
                r_A_1, v_A_1 = fsv.find_orbital_elements(jd_1,self.p_A)[0:2]
                self.r_A = np.append(self.r_A,np.array([r_A_1]))
                r_B_2, v_B_2 = fsv.find_orbital_elements(jd_2,self.p_B)[0:2]
                self.r_B = np.append(self.r_B,np.array([r_B_2]))
                v_dep_A, v_arr_B = PorkchopPlot.lamberts_problem(r_A_1, r_B_2, (jd_2 - jd_1)*24*60*60 , mu_s/(10**9), self.dir, 1)
                self.dvs = np.append(self.dvs,np.linalg.norm(v_dep_A-v_A_1) + np.linalg.norm(v_arr_B-v_B_2))

        # Remove outliers
        self.new_jd_1s = np.array([])
        self.new_jd_2s = np.array([])
        self.new_dvs = np.array([])

        self.min_dv = 100

        for i in range(len(self.jd_1s)):
            if self.dvs[i] < self.min_dv:
                self.min_dv = self.dvs[i]
                self.min_jd_1 = self.jd_1s[i]
                self.min_jd_2 = self.jd_2s[i]
            if self.dvs[i] < self.dv_max:
                self.new_jd_1s = np.append(self.new_jd_1s,self.jd_1s[i])
                self.new_jd_2s = np.append(self.new_jd_2s,self.jd_2s[i])
                self.new_dvs = np.append(self.new_dvs,self.dvs[i])

        # Convert Julian times to DD/MM/YY
        self.depart_dates = [PorkchopPlot.julian_day_number_to_gregorian(int(self.new_jd_1s[i])) for i in range(len(self.new_jd_1s))]
        self.arrive_dates = [PorkchopPlot.julian_day_number_to_gregorian(int(self.new_jd_2s[i])) for i in range(len(self.new_jd_1s))]
        print('### Finished!')

    def get_min_dv(self):
        print(f'Minimum delta-v of {self.min_dv} km/s departing at {PorkchopPlot.julian_day_number_to_gregorian(int(self.min_jd_1))} and arriving at {PorkchopPlot.julian_day_number_to_gregorian(int(self.min_jd_2))}')

    def get_plot(self):
        
        f, (ax1, ax2) = plt.subplots(1, 2, width_ratios=[5, 1])

        ax1.set_title(f'Porkchop plot for {self.p_A}-{self.p_B} Transfer')
        ax1.set_xlabel(f'Departure from {self.p_A} date')
        ax1.set_ylabel(f'Arrival to {self.p_B} date')
        ax1.tick_params(labelsize=5)
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y')) 
        ax1.xaxis.set_major_locator(mdates.DayLocator(interval=40))
        ax1.yaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
        ax1.yaxis.set_major_locator(mdates.DayLocator(interval=60))
        # ax1.set_aspect(aspect=0.3,adjustable='box')
        ax1.scatter(self.depart_dates,self.arrive_dates,c=np.log(self.new_dvs),s=5)
        plt.gcf().autofmt_xdate()
        
        zeros = np.zeros(len(self.new_dvs))
        ax2.set_title('Colour scale')
        ax2.set_ylabel('Delta-v (km/s)')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.set_ticks(np.arange(0, self.dv_max, 2))
        ax2.tick_params(labelsize=5)
        ax2.scatter(zeros,self.new_dvs,c=np.log(self.new_dvs),s=5)

        f.savefig(f'porkchop_{self.p_A}_to_{self.p_B}.png',dpi=1080)
        plt.show()

    def get_dv(self,y1,m1,d1,y2,m2,d2):
        # Get Delta-v of particular departure and arrival date
        departure_get = dt.date(y1,m1,d1)
        arrival_get = dt.date(y2,m2,d2)

        closest_depart = dt.date(1000,1,1)
        for depart in self.depart_dates:
            if abs(departure_get-depart) < abs(departure_get-closest_depart):
                closest_depart = depart

        closest_arrive = dt.date(1000,1,1)
        for i in range(len(self.arrive_dates)):
            if self.depart_dates[i] == closest_depart and abs(arrival_get-self.arrive_dates[i]) < abs(arrival_get-closest_arrive):
                closest_arrive = self.arrive_dates[i]

        print(f'Closest dates: {closest_depart} - {closest_arrive}')

        n = 0
        while self.depart_dates[n] != closest_depart:
            n = n+1
        while self.arrive_dates[n] != closest_arrive:
            n = n+1

        print(f'Delta-v = {self.new_dvs[n]}')

    def plot_3D(self):
        plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlabel(f'{self.p_A} Departure Date')
        ax.set_ylabel(f'{self.p_B} Arrival Date')
        ax.set_zlabel('Delta-v (km/s)')
        ax.scatter(self.jd_1s/(60.0*60*24),self.jd_2s/(60.0*60*24),self.dvs,c=np.log(self.dvs), s = 1)
        plt.savefig(f'porkchop_{self.p_A}_to_{self.p_B}_3d.png')
        plt.show()

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

# Example use of instantiation and instance methods
# Parameters are:           Planet 1, Planet 2,  departure dates range,  #dates to test,  range of flight times (days), #flight times to test, direction,  max dv to plot
earth_to_starman = PorkchopPlot('Earth',   'Starman',  0,1,1,2040, 0,1,1,2042,      600,                 10,24*30,                    50,               "pro",         100)
earth_to_starman.get_plot()

# Example use of class methods
# ut = 0
# d = 1
# m = 7
# y = 2031
# jd = fsv.date_to_JD(ut,d,m,y)
# date = PorkchopPlot.julian_day_number_to_gregorian(int(jd))
# print(f'Compare original date {d}/{m}/{y} to {date.day}/{date.month}/{date.year} date that has been converted to JD and back again ')