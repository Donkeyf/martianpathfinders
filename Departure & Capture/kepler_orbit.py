from asyncio import taskgroups
import numpy as np
import scipy as sp
import orbitalTransforms2 as ot
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import sys
import format_figures as ff
import compute_params as cp

### Jack Naylor's figure formatting setup
ff.startup_plotting(font_size=12)

class KeplerOrbit:
    def __init__(self, name, radius, mu, r_p, i,raan,e,argp,ma,t0,j2):
        self.name = name
        self.radius = radius
        self.mu = mu
        self.r_p = r_p
        self.i = cp.deg_to_rad(i)
        self.raan = cp.deg_to_rad(raan)
        self.e = e
        self.argp = cp.deg_to_rad(argp)
        self.ma = cp.deg_to_rad(ma)
        self.t0 = t0
        self.j2 = j2
        self.a = r_p/(1-e)
        self.n = cp.compute_n(mu,self.a)
        self.period = cp.compute_period(mu,self.a)

        self.ta = None
        self.peri = None
        self.r = None
        self.v_r = None
        self.v_n = None
        self.r_j2 = None
        self.alt = None


    def tle_orbit(self,t,dt):
        # Compute Semi-major Axis
        a = cp.compute_semimajoraxis(self.n,self.mu)

        # Compute the true, eccentric and mean anomalies over this time series
        ta, ea, ma_vec = ot.compute_anomalies(t,self.ma,self.e,self.n)

        # Compute the radial and tangential orbital velocities over this time series
        v_r, v_n = ot.compute_orbital_velocity(a,self.e,ta,self.mu)

        # Convert TLE to the Perifocal Frame
        p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta,a,self.e,self.mu)

        # Create perifocal frame coordinates array
        peri = np.row_stack((p,q,w))

        # Convert Perifocal Frame to ECI Frame
        transform = ot.perifocal_to_eci_matrix(self.i,self.raan,self.argp)    
        # Apply transform to create ECI frame coordinates array
        r = np.dot(transform,peri)

        ### J2 Perturbation Version
        if (dt==None or self.j2==None or self.radius==None):
            self.ta = ta
            self.peri = peri
            self.r = r
            self.v_r = v_r
            self.v_n = v_n
            return ta, a, peri, r, v_r, v_n
        else:

            d_raan = cp.compute_d_raan(self.mu,self.j2,self.radius,self.e,self.a,self.i)
            d_argp = cp.compute_d_argp(self.mu,self.j2,self.radius,self.e,self.a,self.i)

            # print(f'd_raan = {d_raan}, d_argp = {d_argp}')
            raan_vec = np.arange(self.raan,self.raan + d_raan*dt*len(p),d_raan*dt)
            argp_vec = np.arange(self.argp,self.argp + d_argp*dt*len(p),d_argp*dt)
            # print(f'raan {raan_vec[0]} - {raan_vec[-1]} argp {argp_vec[0]} - {argp_vec[0]}')

            x_j2 = np.empty((len(p)))
            y_j2 = np.empty((len(p)))
            z_j2 = np.empty((len(p)))
            dx_j2 = np.empty((len(p)))
            dy_j2 = np.empty((len(p)))
            dz_j2 = np.empty((len(p)))

            for k in range(len(p)):
                x_j2[k],y_j2[k],z_j2[k],dx_j2[k],dy_j2[k],dz_j2[k] = ot.perifocal_to_eci(p[k], q[k], w[k], dp[k], dq[k], dw[k], self.i, raan_vec[k], argp_vec[k])

            
            r_j2 = np.row_stack((x_j2,y_j2,z_j2))
            self.ta = ta
            self.peri = peri
            self.r = r
            self.v_r = v_r
            self.v_n = v_n
            self.r_j2 = r_j2
            self.alt = ot.magnitude(peri) - self.radius

            self.r_ecef = ot.eci_to_ecef(t,r,self.t0)

            return ta, a, peri, r, v_r, v_n, r_j2