import numpy as np
import scipy as sp


def hyperbolic_orbital_elements(r_p, vinf_vec, mu, direction):
    v_inf = np.linalg.norm(vinf_vec)
    
    a = mu/(v_inf ** 2)
    e = r_p/a + 1
    beta = np.arccos(1/e)

    z = [0, 0, 1]

    if direction == 0:
        apse_vec = vinf_vec * np.cos(beta) + np.cross(z, vinf_vec) * np.sin(beta) + z * (np.dot(z, vinf_vec)) * (1 - np.cos(beta))
    else:
        apse_vec = vinf_vec * np.cos(beta) + np.cross(-z, vinf_vec) * np.sin(beta) + -z * (np.dot(-z, vinf_vec)) * (1 - np.cos(beta))

    rp_actual = apse_vec/np.linalg.norm(apse_vec) * r_p

    h_dir = np.cross(rp_actual, v_inf)
    h_unit = h_dir/np.linalg.norm(h_dir)
    h = r_p * np.sqrt(v_inf**2 + (2 * mu)/r_p)
    H = h_unit * h

    N = np.cross([0,0,1], H)

    # Right Ascension of the Ascending Node
    if N[1] >= 0:
        raan = np.arccos(N[0]/n)
    elif N[1] < 0:
        raan = 2*np.pi - np.arccos(N[0]/n)

    i = np.arccos(H[2]/h)

    ta_inf = np.arccos(-1/e)

    if E[2] >= 0:
        argp = np.arccos( np.dot(N,E)/n/e )
    elif E[2] < 0:
        argp = 2*np.pi - np.arccos( np.dot(N,E)/n/e )

    return h, e, a, i, raan, argp, ta_inf
    
    
