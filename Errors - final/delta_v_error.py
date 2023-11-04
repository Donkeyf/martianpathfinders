import numpy as np
import scipy as sp


def delta_v_error(rp_vec, percentage, v_p, mu):
    v_error = np.linalg.norm(v_p) * (1 + percentage)
    r_p = np.linalg.norm(rp_vec)
    v_inf = np.sqrt(v_error ** 2 - (2*mu)/r_p)
    a = mu/(v_inf ** 2)
    e = r_p/a + 1
    beta = np.arccos(1/e)

    z = np.array([0, 0, 1])
    vinf_dir = -(rp_vec * np.cos(beta) + np.cross(z, rp_vec) * np.sin(beta) + z * (np.dot(z, rp_vec - np.cos(beta))))

    vinf_unit = vinf_dir/np.linalg.norm(vinf_dir)

    vinf_vec = vinf_unit * v_inf

    return vinf_vec


v_inf_vec = delta_v_error(np.array([1877803.6559577815, -2509185.6809038236, -1329066.7310026744]), 0.01, np.array([-5273.142169293293, -3230.5932029985156, -1349.98546002923]), 1.327e20)
