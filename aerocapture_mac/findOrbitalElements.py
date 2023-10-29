import numpy as np

def find_orbital_elements(R, V, mu):
    # find the orbital elements when given the state vectors
    # use R1, V1 where applicable for the anomalies at the starting point

    # Distance and Speed
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)
    ## radial velocity
    v_r = np.dot(R,V) / r

    # Specific Angular Momentum
    H = np.cross(R,V)
    h = np.linalg.norm(H)

    # Inclination
    i = np.arccos(H[2]/h)

    # Node Line
    N = np.cross([0,0,1], H)
    n = np.linalg.norm(N)
    
    # Right Ascension of the Ascending Node
    if N[1] >= 0:
        raan = np.arccos(N[0]/n)
    elif N[1] < 0:
        raan = 2*np.pi - np.arccos(N[0]/n)

    # Eccentricity
    E = 1/mu * ((v**2-mu/r)*R - r*v_r*V)
    e = np.linalg.norm(E)
    print(e)

    # Argument of Perigee
    if E[2] >= 0:
        argp = np.arccos( np.dot(N,E)/n/e )
    elif E[2] < 0:
        argp = 2*np.pi - np.arccos( np.dot(N,E)/n/e )
    
    # Semi-major Axis
    a = h**2/mu * 1/(1-e**2)
    

    # Period
    T = 2*np.pi/np.sqrt(mu) * a**(3/2)

    # Mean Motion
    n = 2*np.pi/T

    # Anomalies
    ## True anomaly
    if v_r >= 0:
        ta = np.arccos( np.dot(E/e, R/r) )
    elif v_r < 0:
        ta = 2*np.pi - np.arccos( np.dot(E/e, R/r) )

    ## Eccentric Anomaly
    ea = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(ta/2))

    ## Mean Anomaly
    ma = ea - e*np.sin(ea)

    return h, e, a, T, n, i, raan, argp, ta, ma


