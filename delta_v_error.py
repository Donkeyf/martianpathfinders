import numpy as np
import scipy as sp


def delta_v_error(rp_vec, percentage, v_p, mu):
    v_error = np.linalg.norm(v_p) * (1 + percentage)
    r_p = np.linalg.norm(rp_vec)
    v_inf = np.sqrt(v_error ** 2 - (2*mu)/r_p)
    a = mu/(v_inf ** 2)
    e = r_p/a + 1
    beta = np.arccos(1/e)

    z = [0, 0, 1]
    vinf_dir = -v_p * np.cos(beta) + np.cross(-z, -v_p) * np.sin(beta) + -z * (np.dot(-z, -v_p - np.cos(beta)))

    vinf_unit = vinf_dir/np.linalg.norm(vinf_dir)

    vinf_vec = vinf_unit * v_inf

    return vinf_vec