import numpy as np

# Constants
mu_E = 3.986e5          # Earth gravitational parameter (km^3/s^2)
r_E = 6378              # Earth radius

# Circular parking orbit (500km altitude)
r_Ep = r_E + 750
T_Ep = 2*np.pi/np.sqrt(mu_E) * r_Ep**(3/2)


# Function for delta v required for phase shift
def phasing_deltav(T_0, r, mu, t_error):

    # time difference
    t_AB = t_error * 86400

    # calculate orbital parameters for phasing orbit
    r_p_1 = r
    T_1 = T_0 - t_AB
    a_1 = (T_1*np.sqrt(mu)/2/np.pi)**(2/3)
    r_a_1 = 2*a_1 - r_p_1
    h_1 = np.sqrt(2*mu)*np.sqrt(r_a_1*r_p_1/(r_a_1+r_p_1))

    # calculate required delta v
    v_0 = np.sqrt(mu/r)
    v_1 = h_1/r_p_1
    dv = 2*abs(v_0-v_1)

    return dv


# Time error
t_error1 = 0.01 / 100
t_error2 = 0.05 / 100
t_error3 = 0.1 / 100


dv1 = phasing_deltav(T_Ep, r_Ep, mu_E, t_error1)
dv2 = phasing_deltav(T_Ep, r_Ep, mu_E, t_error2)
dv3 = phasing_deltav(T_Ep, r_Ep, mu_E, t_error3)





print(dv1)
print(dv2)
print(dv3)