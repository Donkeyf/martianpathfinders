import numpy as np
import scipy as sp

def find_state_vector(v_inf, r_p, mu, soi):
    v_mag = np.linalg.norm(v_inf)
    a = mu/(v_inf**2)
    e = r_p/a + 1
    delta = a * np.sqrt(e**2 - 1)

    asymptote_vec = -v_inf/v_mag * soi

    