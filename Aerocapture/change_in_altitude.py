import numpy as np
import scipy as sp
import orbitalTransforms as ot
import computeHyperbolicAnomalies as cha
import matplotlib.pyplot as plt

mu_mars = 4.282831567e13
r_mars = 3396.2e3


def altitude(v_inf, mars_radius, mu, t_vec):
    a_h = mu/(v_inf**2)
    e_h = ((mars_radius + 80000)/a_h) + 1
    h_h = np.sqrt(a_h * mu * ((e_h**2) - 1))
    ta_h, ea, ma_vec = cha.compute_hyperbola_anomalies(t_vec, -np.pi/2, e_h, h_h, mu)

    p_h, q_h, w_h, dp_h, dq_h, dw_h = ot.elements_to_perifocal(ta_h, a_h, e_h, mu, h_h)
    x_h, y_h, z_h, dx_h, dy_h, dz_h = ot.perifocal_to_eci(p_h, q_h, w_h, dp_h, dq_h, dw_h, 0, 0, 0)
    
    altitude = []
    i = 0
    min = 1000000000
    min_step = 0
    while i < len(x_h):
        temp = np.linalg.norm([x_h[i], y_h[i], z_h[i]]) - 3396.2e3
        altitude.append(temp)
        if temp < min:
            min = temp
            min_step = i
        i += 1
    
    print(min)
    print(min_step)

    return altitude

t_vec = np.linspace(0, 600, 600)
alt = altitude(6383, r_mars, mu_mars, t_vec)

print(alt[-1])

plt.figure(1)
plt.plot(t_vec, alt)
plt.show()