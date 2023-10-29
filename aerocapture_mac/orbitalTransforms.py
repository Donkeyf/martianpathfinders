import numpy as np
import scipy as sp


def compute_anomalies(t_sec, ma, e, n=0):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity

    # Propagate mean anomaly using mean motion (vector over time)
    ma_vec = t_sec * n + ma

    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(E):
        return ma_vec - E + float(e) * np.sin(E)

    init_guess = ma_vec

    # Calc eccentric anomaly
    ea = sp.optimize.newton(kepler_equation, init_guess)

    # Calc true anomaly
    ta = 2 * np.arctan(((1 + e)/(1 - e))**(1/2) * np.tan(ea/2))

    # Wrap anomalies to -pi:pi
    ta = [i - 2* np.pi if i>np.pi else i for i in ta]
    ea = [i - 2* np.pi if i>np.pi else i for i in ea]
    ma_vec = [i - 2* np.pi if i>np.pi else i for i in ma_vec]
    return ta, ea, ma_vec


def compute_orbital_velocity(h, e, ta, mu):

    v_r = mu/h * e * np.sin(ta)   # Radial velocity
    v_n = mu/h * (1 + e * np.cos(ta)) # Normal/Tangential velocity

    return v_r, v_n


def elements_to_perifocal(ta, a, e, mu, h):
    # Calc perifocal distance in m
    r = (h ** 2)/mu * 1/(1 + e * np.cos(ta))

    # Compute perifocal coordinates
    p = r * np.cos(ta)
    q = r * np.sin(ta)
    w = 0

    # Compute perifocal velocities
    dp = mu/h * -np.sin(ta)
    dq = mu/h * (e + np.cos(ta))
    dw = 0

    return p, q, w, dp, dq, dw


def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):


    # Transform coordinates to ECI frame
    
    x = np.zeros_like(p)
    y = np.zeros_like(p)
    z = np.zeros_like(p)

    dx = np.zeros_like(p)
    dy = np.zeros_like(p)
    dz = np.zeros_like(p)

    # Compute coordinates
    k = 0
    print(len(x))
    while k < len(x):
        r_vec = np.matmul(perifocal_to_eci_matrix(i, raan, argp), [[p[k]], [q[k]], [w]])

        x[k] = r_vec[0]
        y[k] = r_vec[1]
        z[k] = r_vec[2]
        k += 1

    j = 0
    while j < len(dx):
        v_vec = np.matmul(perifocal_to_eci_matrix(i, raan, argp), [[dp[j]], [dq[j]], [dw]])
        dx[j] = v_vec[0]
        dy[j] = v_vec[1]
        dz[j] = v_vec[2]
        j += 1

    return x, y, z, dx, dy, dz


def perifocal_to_eci_matrix(i, raan, argp):

    # Calculate transformation matrix from perifocal to ECI frame
    

    p_to_e = np.matrix([[np.cos(raan)*np.cos(argp)-np.sin(raan)*np.cos(i)*np.sin(argp), -np.cos(raan)*np.sin(argp)-np.sin(raan)*np.cos(i)*np.cos(argp), np.sin(raan)*np.sin(i)],
             [np.sin(raan)*np.cos(argp)-np.cos(raan)*np.cos(i)*np.sin(argp), -np.sin(raan)*np.sin(argp)+np.cos(raan)*np.cos(i)*np.cos(argp), np.cos(raan)*np.sin(i)], 
             [np.sin(i)*np.sin(argp), np.sin(i)*np.cos(argp), np.cos(i)] 
              ])

    return p_to_e
