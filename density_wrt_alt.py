import numpy as np
import scipy as sp

#Atmospheric models for aerocapture

def mars_atmosphere_model(alt):
    R = 191.181
    r0 = 3389.51
    h = alt * 0.001
    z = (-h * r0)/(h - r0)
    


    if  h < 39:
        P = 610.5 * (228.5/(228.50 - 1.8 * h)) ** (-19.435/1.8)
        T = 226.5 - 1.8 * h
    elif 39 < h < 48:
        P = 11.6025 * np.exp(-19.35 * (h - 39) / 158.30)
        T = 158.30
    elif 48 <= h < 55:
        T = 271.10 - 2.35 * h
        P = 3.84305 * (158.30 / (158.30 - 2.35 * (h - 48))) ** (19.435/-2.35)
    elif 55 <= h < 66:
        T = 106.10 + 0.65 * h
        P = 1.55091 * (141.85 / (141.85 + 0.65 * (h - 55))) ** (19.435/0.65)
    elif 66 <= h < 75:
        T = 314.00 - 2.50 * h
        P = 0.356464 * (149.00 / (149.00 - 2.50 * (h - 66))) ** (19.435 / -2.50)
    elif 75 <= h < 84:
        T = -61.00 + 2.50 * h
        P = 0.0998430 * (126.50 / (126.50 + 2.50 * (h - 75))) ** (19.435 / 2.50)
    elif 84 <= h < 95:
        T = 149.00
        P = 0.0279653 * np.exp(-19.435 * (h - 84) / 149.00)
    elif 95 <= h < 105:
        T = 282.00 - 1.40 * h
        P = 0.00666032 * (149.00 / (149.00 - 1.40 * (h - 95))) ** (19.435 / -1.40)
    elif 105 <= h < 115:  
        T = 203.25 - 0.65 * h
        P = 0.00169282 * (135.00 / (135.00 - 0.65 * (h - 105))) ** (-19.435/0.65)
    elif 115 <= h < 200:
        temp = (-2.55314E-10 * z ** 5 + 2.31927E-07 * z **4 - 8.33206E-05 * z ** 3 + 0.0151947 * z ** 2 - 1.52799 * z + 48.69659)
        return np.exp(temp)
    elif 200 <= h < 300:
        temp = (2.65472E-11 * z ** 5 - 2.45558E-08 * z **4 + 6.31410E-06 * z ** 3 + 4.73359E-04 * z ** 2 - 0.443712 * z + 23.79408)
        return np.exp(temp)
    else:
        return 0

    
    rho = P/(R * T)

    return rho

def earth_atmosphere_model(h):
    # R = 287.053
    # h = h * 0.001

    # if 20 < h < 32:
    #     T = 196.65 + h
    #     P = 5474.889 * (216.65 / (216.65 + (h - 20))) ** (34.1632)
    # elif 32 <= h < 47:
    #     T = 139.05 + 2.8 * h
    #     P = 868.0187 * (228.65 / (228.65 + 2.8 * (h - 32))) ** (34.1632 / 2.8)
    # elif 47 <= h < 51:
    #     T = 270.65
    #     P = 110.9063 * np.exp(-34.1632 * (h - 47) / 270.65)
    # elif 51 <= h < 71:
    #     T = 413.45 - 2.8 * h  
    #     P = 66.93887 * (270.65 / (270.65 - 2.8 * (h - 51))) ** (34.1632 / -2.8)  
    # elif 71 <= h < 85:
    #     T = 356.65 - 2.0 * h
    #     P = 3.956420 * (214.65 / (214.65 - 2 * (h - 71))) ** (34.1632 / -2)
    # elif 85 <= h < 91:
    #     b = -3.322622E-06
    #     c = 9.111460E-04
    #     d = -0.2609971
    #     e = 5.944694
    # elif 91 <= h < 100:
    #     b = 2.873405E-05
    #     c = -0.008492037
    #     d = 0.6541179
    #     e = -23.62010
    # elif 100 <= h < 110:
    #     b = 0.005162063
    #     c = -0.8048342
    #     d = 55.55996
    #     e = -1443.338
    # elif 110 <= h < 120:
    #     b = -8.854164E-05
    #     c = 0.03373254
    #     d = -4.390837
    #     e = 176.5294

    # if 32 < h < 85:
    #     return P/(R * T)
    # elif 85 <= h < 120:
    #     r0 = 6356.766 - 6378
    #     z = (-h * r0)/(h - r0)
    #     # print(np.exp(b * z ** 3 + c * z ** 2 + d * z + e))
    #     return np.exp(b * z ** 3 + c * z ** 2 + d * z + e)
    # else:
    #     return 0

    if h <25000:
        T = -56.46
        p = 22.65 * np.exp(1.73-0.000157 * h)
    elif 25000 < h:
        T = -131.21 + 0.00299 * h
        p = 2.488 * ((T + 273.1)/216.6) ** -11.388
    else:
        return 0
    
    return p/(0.28969 * (T + 273.1))