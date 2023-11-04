import numpy as np
import scipy as sp
import density_wrt_alt as d_alt


#aerobrake numerical solver function for mars
#y0 - initial value ([x, y, z, dx, dy, dz])
#tspan - time span ([t0, tfinal])
#t_eval - solution time steps (has to be within tspan)
#m - spacecraft mass
#orbital_body - 0 = mars, 1 = earth
def aero_capture(y0, t_span, t_eval, m, mu, radius, orbital_body):

    def drag(t, f):
        r_vec = f[:3] #position
        r_dot = f[3:] #velocity

        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(r_dot)

        h = r - radius
        # T = -23.4 - 0.00222 * h
        # p = 0.699 * np.exp(-0.00009 * h) 
        # rho = p/(0.1921 * (T + 273.1))
        if orbital_body == 0:    #mars
            rho = d_alt.mars_atmosphere_model(h)
            r_ddot = -mu/(r**3) * r_vec - (11.8124 * rho * v * r_dot)/m
        elif orbital_body == 1:  #earth
            rho = d_alt.earth_atmosphere_model(h)
            r_ddot = -mu/(r**3) * r_vec - (0.52875 * rho * v * r_dot)/m

        #two body equation with drag term
        r_ddot = -mu/(r**3) * r_vec - (11.8124 * rho * v * r_dot)/m

        
        return np.concatenate((r_dot, r_ddot))


    # return sp.integrate.RK45(drag, 0, y0, 2223, rtol=1e-6, atol=1e-6, vectorized=True)
    if orbital_body == 0:
        return sp.integrate.solve_ivp(drag, t_span, y0, t_eval=t_eval, dense_output=True, rtol=1e-5, atol=1e-5)
    else:
        return sp.integrate.solve_ivp(drag, t_span, y0, t_eval=t_eval, dense_output=True, rtol=1e-4, atol=1e-4)
