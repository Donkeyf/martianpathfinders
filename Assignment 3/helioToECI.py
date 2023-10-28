import numpy as np

def helio_to_eci(Vec_Helio, planet="Earth"):

    # convert vectors between Heliocentric frame to ECI frame
    # only applies for velocities, not positions 

    # obliquity
    if planet == "Earth":
        o = 23.44 * np.pi/180
    elif planet == "Mars":
        o = 25.19 * np.pi/180

    # matrix rotation about the x-axis
    R_mat = np.matrix([[1, 0, 0], [0, np.cos(o), np.sin(o)], [0, -np.sin(o), np.cos(o)]])
    Vec_ECI = np.array(np.matmul(R_mat, Vec_Helio))[0]

    # calculate inclination angle of excess velocity
    i_calc = abs( np.arcsin(Vec_ECI[2]/np.linalg.norm([Vec_ECI[0],Vec_ECI[1],Vec_ECI[2]])) )

    return Vec_ECI, i_calc


# Miscellaneous: rotation matrices about other axes
phi = 0
R_mat_y = np.matrix([[np.cos(phi), 0, -np.sin(phi)], [0, 1, 0], [np.sin(phi), 0, np.cos(phi)]])
R_mat_z = np.matrix([[np.cos(phi), np.sin(phi), 0], [-np.sin(phi), np.cos(phi), 0], [0, 0, 1]])