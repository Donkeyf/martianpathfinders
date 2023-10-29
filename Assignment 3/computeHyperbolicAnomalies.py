import numpy as np

def compute_hyperbola_anomalies(t_vec, ma, e, h, mu):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity

    # Propagate mean anomaly
    ma_vec = t_vec * mu**2/h**3 * (e**2-1)**(3/2) + ma

    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(F):
        return e * np.sinh(F) - F

    # Calculate eccentric anomaly --- Newton's Method
    ea = np.zeros(len(ma_vec))      # initialise eccentric anomaly vector

    for j in range(len(ma_vec)):
        F0 = ma_vec[j]              # set initial guess for E
        tol = 1                     # initialise tolerance value
   
        while tol > 1e-8:           # iterate until tolerance condition (<1e-8) is satisfied
            f = kepler_equation(F0) - ma_vec[j]
            f_prime = e * np.cosh(F0) - 1
            F1 = F0 - f/f_prime     # Newton's method

            tol = abs(F1-F0)        # calculate tolerance
            F0 = F1                 # update guess value for following iteration

        ea[j] = F1                  # update eccentric anomaly vector

    # Calculate true anomaly
    ta = 2 * np.arctan( np.sqrt((e+1)/(e-1)) * np.tanh(ea/2) )

    # Wrap anomalies to -pi:pi
    ta = [i-2*np.pi if i>np.pi else i for i in ta]
    ea = [i-2*np.pi if i>np.pi else i for i in ea]
    ma_vec = [i-2*np.pi if i>np.pi else i for i in ma_vec]

    return ta, ea, ma_vec