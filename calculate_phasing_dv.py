import numpy as np

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

