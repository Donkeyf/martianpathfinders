import numpy as np
import scipy as sp
import orbitalTransforms2 as ot
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import sys
import format_figures as ff

def soi(R,mp,ms):
    return R* ((mp/ms)**(2/5))

def deg_to_rad(angle):
    angle = angle/180.0 * np.pi
    return angle

def compute_semimajoraxis(n,mu):
    return np.cbrt(mu/(n**2))

def convert_n(n_day):
    return (n_day*2*np.pi)/(60*60*24)

def convert_ta_to_ma(ta,e):
    ea = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(ta/2))
    return ea - e*np.sin(ea)

def compute_h(mu,ra,rp):
    return np.sqrt(2*mu)*np.sqrt((ra*rp)/(ra+rp))

def compute_e(ra,rp):
    return (ra-rp)/(ra+rp)

def compute_v(mu,r,a):
    return np.sqrt(2*mu)*np.sqrt(1/r - 1/(2*a))

def compute_period(mu,a):
    return (2*np.pi)/np.sqrt(mu) * np.sqrt(a**3)

def compute_n(mu,a):
    return np.sqrt(mu/a**3)

def set_aspect_equal(ax,x,y,z):
    min = np.min(np.array([np.min(x),np.min(y),np.min(z)]))
    max = np.max(np.array([np.max(x),np.max(y),np.max(z)]))
    ax.set_xlim3d(min,max)
    ax.set_ylim3d(min,max)
    ax.set_zlim3d(min,max)

def compute_an_dn(i,raan,e,mu,a,argp):
    # AN is at ta = -argp
    p_an, q_an, w_an, dp_an, dq_an, dw_an = ot.elements_to_perifocal([-argp],a,e,mu)
    peri_an = np.row_stack((p_an,q_an,w_an))

    ### Convert Perifocal Frame to ECI Frame
    transform = ot.perifocal_to_eci_matrix(i,raan,argp)    
    # Apply transform to create ECI frame coordinates array
    r_an = np.dot(transform,peri_an)

    # DN is at ta = pi-argp
    p_dn, q_dn, w_dn, dp_dn, dq_dn, dw_dn = ot.elements_to_perifocal([np.pi-argp],a,e,mu)
    peri_dn = np.row_stack((p_dn,q_dn,w_dn))

    # Apply transform to create ECI frame coordinates array
    r_dn = np.dot(transform,peri_dn)

    return p_an, q_an, w_an, r_an,  p_dn, q_dn, w_dn, r_dn

"""
dt, j2, r_p - optional parameters for J2 perturbation
"""
def tle_orbit(t, mu, i,raan,e,argp,ma,n,dt="not given",j2="not given",r_p="not given"):
    # Compute Semi-major Axis
    a = compute_semimajoraxis(n,mu)

    # Compute the true, eccentric and mean anomalies over this time series
    ta, ea, ma_vec = ot.compute_anomalies(t,ma,e,n)

    # Compute the radial and tangential orbital velocities over this time series
    v_r, v_n = ot.compute_orbital_velocity(a,e,ta,mu)

    # Convert TLE to the Perifocal Frame
    p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta,a,e,mu)

    # Create perifocal frame coordinates array
    peri = np.row_stack((p,q,w))

    # Convert Perifocal Frame to ECI Frame
    transform = ot.perifocal_to_eci_matrix(i,raan,argp)    
    # Apply transform to create ECI frame coordinates array
    r = np.dot(transform,peri)

    ### J2 Perturbation Version
    if (dt=="not given" or j2=="not given" or r_p=="not given"):
        return ta, a, peri, r, v_r, v_n
    else:
        d_raan = compute_d_raan(mu,j2,r_p,e,a,i)
        d_argp = compute_d_argp(mu,j2,r_p,e,a,i)

        print(f'd_raan = {d_raan}, d_argp = {d_argp}')

        raan_vec = np.arange(raan,raan + d_raan*dt*len(p),d_raan*dt)
        argp_vec = np.arange(argp,argp + d_argp*dt*len(p),d_argp*dt)
        print(f'raan {raan_vec[0]} - {raan_vec[-1]} argp {argp_vec[0]} - {argp_vec[0]}')

        x_j2 = np.empty((len(p)))
        y_j2 = np.empty((len(p)))
        z_j2 = np.empty((len(p)))
        dx_j2 = np.empty((len(p)))
        dy_j2 = np.empty((len(p)))
        dz_j2 = np.empty((len(p)))

        for k in range(len(p)):
            x_j2[k],y_j2[k],z_j2[k],dx_j2[k],dy_j2[k],dz_j2[k] = ot.perifocal_to_eci(p[k], q[k], w[k], dp[k], dq[k], dw[k], i, raan_vec[k], argp_vec[k])

        
        r_j2 = np.row_stack((x_j2,y_j2,z_j2))
        return ta, a, peri, r, v_r, v_n, r_j2


def compute_d_raan(mu,j2,r_p,e,a,i):
    return -1 * ((3*np.sqrt(mu)*j2*r_p**2)/(2*(1-e**2)**2 * a**(7/2))) * np.cos(i)

def compute_d_argp(mu,j2,r_p,e,a,i):
    return -1 * ((3*np.sqrt(mu)*j2*r_p**2)/(2*(1-e**2)**2 * a**(7/2))) * (5/2 * np.sin(i)**2 - 2)