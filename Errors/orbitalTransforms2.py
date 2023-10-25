import numpy as np
import scipy as sp


def compute_anomalies(t_sec, ma, e, n=0):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity

    # Propagate mean anomaly using mean motion (vector over time)
    ma_vec = ma + (n * t_sec) # mean anomaly at time t_sec is the initial mean anomaly ma plus the mean motion n times t_sec

    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(E):
        return ma_vec - (E - e * np.sin(E)) # this gives mean anomaly as defined in the lectures

    init_guess = ma_vec

    # Calc eccentric anomaly
    ea = sp.optimize.newton(kepler_equation,init_guess)

    # Calc true anomaly
    
    ta = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(ea/2))

    def wrap(angles):
        angles = angles%(2*np.pi)
        return angles*(1-(angles//np.pi)) + (angles-2*np.pi)*((angles//np.pi))
        # return angles
    # angles%(2*np.pi) 0 to 2pi
    # (angles//np.pi)  0 = 0 to pi, 1 = pi to 2pi       MASK for pi to 2pi
    # 1-(angles//np.pi) 1 = 0 to pi, 0 = pi to 2pi      MASK for 0 to pi

    # Wrap anomalies to -pi:pi
    ta = wrap(ta)
    ea = wrap(ea)
    ma_vec = wrap(ma_vec)
    return ta, ea, ma_vec


def compute_orbital_velocity(a, e, ta, mu):

    h = np.sqrt(a*mu*(1-e**2))  # Specific angular momentum

    v_r = (mu/h) * e*np.sin(ta)  # Radial velocity
    v_n = (mu/h) * (1 + e*np.cos(ta)) # Normal/Tangential velocity

    return v_r, v_n

def elements_to_perifocal(ta, a, e, mu, h="not given"):

    if h=="not given":
        h = np.sqrt(a*mu*(1-e**2))  # Specific angular momentum
    
    # Calc perifocal distance in m
    r = (h**2/mu) * 1/(1+e*np.cos(ta))

    # Compute perifocal coordinates
    p = r*np.cos(ta)
    q = r*np.sin(ta)
    w = np.zeros(len(ta))

    # Compute perifocal velocities
    dp = (mu/h) * (-1) * np.sin(ta)
    dq = (mu/h) * (e + np.cos(ta))
    dw = np.zeros(len(ta))

    return p, q, w, dp, dq, dw


def magnitude(vector):
    return np.sqrt(sum(element**2 for element in vector))

def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):
    # Transform coordinates to ECI frame

    an = np.array([np.cos(raan),np.sin(raan),0]) # Ascending Node unit vector
    zu = np.array([0,0,1]) # z-axis unit vector

    wu = np.cos(i)*zu + np.sin(i)*np.cross(an,zu)
    pu = np.cos(argp)*an + np.sin(argp)*np.cross(wu,an)
    qu = np.cross(wu,pu)

    # print(f'an is {an}, length = {magnitude(an)}')
    # print(f'zu is {zu}, length = {magnitude(zu)}')
    # print(f'wu is {wu}, length = {magnitude(wu)}')
    # print(f'pu is {pu}, length = {magnitude(pu)}')
    # print(f'qu is {qu}, length = {magnitude(qu)}')

    x = np.zeros_like(p)
    y = np.zeros_like(p)
    z = np.zeros_like(p)

    dx = np.zeros_like(p)
    dy = np.zeros_like(p)
    dz = np.zeros_like(p)

    # Compute coordinates
    r = p*pu + q*qu + w*wu
    x = r[0]
    y = r[1]
    z = r[2]


    v = dp*pu + dq*qu + dw*wu
    dx = v[0]
    dy = v[1]
    dz = v[2]

    # Matrix method
    # transform = perifocal_to_eci_matrix(i,raan,argp)

    # peri = np.column_stack((p,q,w))
    # eci = np.zeros_like(peri)
    # d_peri = np.column_stack((dp,dq,dw))
    # d_eci = np.zeros_like(d_peri)
    # for i in range(len(peri)):
    #     eci[i] = np.dot(transform,peri[i])
    #     x[i] = eci[i][0]
    #     y[i] = eci[i][1]
    #     z[i] = eci[i][2]

    #     dx[i] = d_eci[i][0]
    #     dy[i] = d_eci[i][1]
    #     dz[i] = d_eci[i][2]
    

    return x, y, z, dx, dy, dz


def perifocal_to_eci_matrix(i, raan, argp):

    # Calculate transformation matrix from perifocal to ECI frame
    
    an = np.array([np.cos(raan),np.sin(raan),0]) # Ascending Node unit vector
    zu = np.array([0,0,1]) # z-axis unit vector

    wu = np.cos(i)*zu + np.sin(i)*np.cross(an,zu)
    pu = np.cos(argp)*an + np.sin(argp)*np.cross(wu,an)
    qu = np.cross(wu,pu)

    p_to_e = np.linalg.inv(np.row_stack((pu,qu,wu)))

    return p_to_e

def eci_to_ecef(t,r,t0):
    
    # Convert TLE time to Y/M/D
    date = np.datetime64("20" + t0[0:2],'Y') + np.timedelta64(int(float(t0[2:])*24*60*60*1000),'ms')
    # print(date)
    y = date.astype('datetime64[Y]').astype(int) + 1970
    m = date.astype('datetime64[M]').astype(int) % 12 + 1
    d = (date.astype('datetime64[D]') - date.astype('datetime64[M]').astype('datetime64[D]')).astype(int) + 1
    # print(f'{y}/{m}/{d}')

    # Convert epoch to Julian days
    j0 = (367*y - int((7*(y+int((m+9)/12)))/4) + int(275*m/9) + d) + 1721013.5
    # print(j0)

    # Convert epoch
    t0 = (j0-2451545)/36525

    # Convert to degrees
    thetaG0 = 100.4606184 + 36000.77004*t0 + 0.000387933*t0**2 - 2.583*(10**(-8))*t0**3
    
    # Calculate universal time
    hours = (date.astype('datetime64[h]') - date.astype('datetime64[D]').astype('datetime64[h]')).astype(int)
    minutes = (date.astype('datetime64[m]') - date.astype('datetime64[h]').astype('datetime64[m]')).astype(int)
    seconds = (date.astype('datetime64[s]') - date.astype('datetime64[m]').astype('datetime64[s]')).astype(int)
    # print(f'{hours}:{minutes}:{seconds}')
    ut = hours + minutes/60 + seconds/3600
    # print(ut)

    theta_G = ((thetaG0 + 360.98564724 * ut/24) % 360)/180.0 * np.pi
    # print(theta_G)

    w_e = 7.2921150 * 10**(-5) # angular speed of the Earth in inertial space rad/s
    # t_init = 11402014  # start epoch - vernal equinox (OUTDATED: USING POLYNOMIAL NOW)
    theta = theta_G + w_e*t
    r = np.transpose(r)
    r_ecef = np.empty(r.shape)

    for i in range(len(r)):
        rot_mat = np.array([[np.cos(theta[i]),np.sin(theta[i]),0] , [-1*np.sin(theta[i]),np.cos(theta[i]),0] , [0,0,1]])
        r_ecef[i] = np.dot(rot_mat,r[i])

    r_ecef = np.transpose(r_ecef)

    return r_ecef



