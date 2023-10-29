import numpy as np
import scipy as sp


#r_p - radius of perigee
#vinf_vector - vector of the excess velocity
#mu - standard gravitational parameter
#direction:
# 0 =  arrival from soi
# 1 = departure from soi

def hyperbolic_orbital_elements(r_p, vinf_vec, mu, direction):
    v_inf = np.linalg.norm(vinf_vec)
    print(v_inf)
    vinf_vec = np.array(vinf_vec)

    
    a = mu/(v_inf ** 2)
    e = r_p/a + 1
    beta = np.arccos(1/e)

    z = np.array([0, 0, 1])

    if direction == 0:   #Arrival: direction = 0
        apse_vec = (vinf_vec * np.cos(beta) + np.cross(z, vinf_vec) * np.sin(beta) + z * (np.dot(z, vinf_vec)) * (1 - np.cos(beta)))
    else:                #Departure: direction = 1
        apse_vec = -(vinf_vec * np.cos(beta) + np.cross(z, vinf_vec) * np.sin(beta) + z * (np.dot(z, vinf_vec)) * (1 - np.cos(beta)))

    #periapse position vector
    rp_actual = apse_vec/np.linalg.norm(apse_vec) * r_p

    #specific angular momentum
    h_dir = np.cross(rp_actual, vinf_vec)
    h_unit = h_dir/np.linalg.norm(h_dir)
    h = r_p * np.sqrt(v_inf**2 + (2 * mu)/r_p)
    H = h_unit * h

    #node line
    N = np.cross([0,0,1], H)

    # Right Ascension of the Ascending Node
    n = np.linalg.norm(N)
    if N[1] >= 0:
        raan = np.arccos(N[0]/n)
    elif N[1] < 0:
        raan = 2*np.pi - np.arccos(N[0]/n)

    #inclination
    i = np.arccos(H[2]/h)

    #ta at soi
    ta_inf = np.arccos(-1/e)

    #eccentricity
    E = rp_actual/np.linalg.norm(rp_actual) * e

    #argument of perigee
    if E[2] >= 0:
        argp = np.arccos( np.dot(N,E)/n/e )
    elif E[2] < 0:
        argp = 2*np.pi - np.arccos( np.dot(N,E)/n/e )

    return h, e, a, i, raan, argp, ta_inf
    
    
