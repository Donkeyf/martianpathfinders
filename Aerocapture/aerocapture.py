import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import orbitalTransforms as ot
import findOrbitalElements as oe


#aerocapture transfer orbit
mu_mars = 4.282831567e13
r_mars = 3396.2e3
ra = r_mars + 400000
rp = r_mars + 500000
h_t = np.sqrt(2*mu_mars)*np.sqrt(ra*rp/(ra + rp))
v_p = h_t/rp



#aerobrake numerical solver
def aero_capture(y0, t_span):
    def drag(t, f):
        mu_mars = 4.282831567e13
        r_mars = 3396.2e3
        m = 31700
        r_vec = f[:3] #position
        r_dot = f[3:6] #velocity

        # print(f[6])
        r = np.linalg.norm(r_vec)
        #print(r)
        v = np.linalg.norm(r_dot)

        r_ddot = -mu_mars/(r**3) * r_vec - (1.6611 * f[6] * v * r_dot)/m
        h = r - r_mars
        
        # T = -23.4 - 0.00222 * h
        # p = 0.699 * np.exp(-0.00009 * h) 
        # rho_dot = p/(0.1921 * (T + 273.1))

        add = np.concatenate((r_dot, r_ddot))
        test = np.concatenate((add, [rho_dot]))
        #print(test)
        return test

    # return sp.integrate.RK45(drag, 0, y0, 2223, rtol=1e-6, atol=1e-6, vectorized=True)
    return sp.integrate.solve_ivp(drag, t_span, y0, dense_output=True, rtol=1e-6, atol=1e-6)


def hyperbolic_trajectory_state_vector(v_inf, alt, mars_radius, mu):
    a_h = mu/(v_inf**2)
    e_h = ((mars_radius + 80000)/a_h) + 1
    print(e_h)
    h_h = np.sqrt(a_h * mu * ((e_h**2) - 1))
    theta_inf = np.arccos(-1/e_h)
    ta_h = np.arange(-theta_inf, 0, np.pi/180)
    p_h, q_h, w_h, dp_h, dq_h, dw_h = ot.elements_to_perifocal(ta_h, a_h, e_h, mu, h_h)
    x_h, y_h, z_h, dx_h, dy_h, dz_h = ot.perifocal_to_eci(p_h, q_h, w_h, dp_h, dq_h, dw_h, 0, 0, 0)

    
    return [x_h[50], y_h[50], z_h[50], dx_h[50], dy_h[50], dz_h[50]]
    i = 0
    while True:
        l = np.linalg.norm([x_h[i], y_h[i], z_h[i]])
        if np.linalg.norm([x_h[i], y_h[i], z_h[i]]) < alt:
            return [x_h[i], y_h[i], z_h[i], dx_h[i], dy_h[i], dz_h[i]]      
        i += 1


y0 = hyperbolic_trajectory_state_vector(6383.46, 130000 + r_mars, r_mars, mu_mars)
state = np.concatenate((y0, [0]))

t_span = [0, 10000]
dragpath = aero_capture(state, t_span)

#h, e, a, T, n, i, raan, argp, ta, ma = oe.find_orbital_elements([dragpath.y[0], dragpath.y[1], dragpath.y[2]], )


#Plotting Mars
r_mars = 3396.2e3
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_mars = r_mars * np.outer(np.cos(u), np.sin(v))
y_mars = r_mars * np.outer(np.sin(u), np.sin(v))
z_mars = r_mars * np.outer(np.ones(np.size(u)), np.cos(v))


plt.figure(4)
ax = plt.axes(projection='3d')
ax.plot3D(dragpath.y[0], dragpath.y[1], dragpath.y[2])
#ax.plot3D(y0[0][10:], y0[1][10:], y0[2][10:])
ax.plot_surface(x_mars, y_mars, z_mars, color='r')
ax.set_title('')
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_aspect('equal')
plt.show()

# kms = []
# i = 0
# while i < len(dragpath.y[3]):
#     kms.append(np.linalg.norm([dragpath.y[3][i], dragpath.y[4][i], dragpath.y[5][i]]))
#     i += 1

# t = np.arange(0, 190, 1)
# plt.figure(5)
# plt.plot(t, kms)
# plt.show()