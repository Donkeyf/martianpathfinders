import numpy as np
import scipy as sp

#Atmospheric models for aerocapture

def mars_atmosphere_model(h):
    R = 191.181
    r0 = 3389.51
    z = (-h * r0)/(h - r0)

    if  h < 39000:
        P = 610.5 * (228.5/(228.50 - 1.8 * h)) ** (19.435/-1.8)
        T = 226.5 - 1.8 * h
    elif 39000 < h < 48000:
        11.6025 * np.exp(-19.35 * (h - 39) / 158.30)
        T = 158.30
    elif 48000 <= h < 55000:
        T = 271.10 - 2.35 * h
        P = 3.84305 * (158.30 / (158.30 - 2.35 * (h - 48))) ** (19.435/-2.35)
    elif 55000 <= h < 66000:
        T = 106.10 + 0.65 * h
        P = 1.55091 * (141.85 / (141.85 + 0.65 * (h - 55))) ** (19.435/0.65)
    elif 66000 <= h < 75000:
        T = 314.00 - 2.50 * h
        P = 0.356464 * (149.00 / (149.00 - 2.50 * (h - 66))) ** (19.435 / -2.50)
    elif 75000 <= h < 84000:
        T -61.00 + 2.50 * h
        P = 0.0998430 * (126.50 / (126.50 + 2.50 * (h - 75))) ** (19.435 / 2.50)
    elif 84000 <= h < 95000:
        T = 149.00
        P = 0.0279653 * np.exp(-19.435 * (h - 84) / 149.00)
    elif 95000 <= h < 105000:
        T = 282.00 - 1.40 * h
        P = 0.00666032 * (149.00 / (149.00 - 1.40 * (h - 95))) ** (19.435 / -1.40)
    elif 105000 <= h < 115000:  
        T = 203.25 - 0.65 * h
        P = 0.00169282 * (135.00 / (135.00 - 0.65 * (h - 105))) ** (19.435 / -0.65)
    elif 115000 <= h < 200000:
        return np.exp(-2.55314E-10 * z ** 5 + 2.31927E-07 * z **4 - 8.33206E-05 * z ** 3 + 0.0151947 * z ** 2 - 1.52799 * z + 48.69659)
    elif 200000 <= h < 300000:
        return np.exp(2.65472E-11 * z ** 5 - 2.45558E-08 * z **4 + 6.31410E-06 * z ** 3 + 4.73359E-04 * z ** 2 - 0.443712 * z + 23.79408)
    
    else:
        P = 0
        T = 1
    
    rho = P/(R * T)

    return rho

def earth_atmosphere_model(h):
    R = 287.053

    if 32000 < h < 47000:
        T = 139.05 + 2.8 * h
        P = 868.0187 * (228.65 / (228.65 + 2.8 * (h - 32))) ** (34.1632 / 2.8)
    elif 47000 <= h < 51000:
        T = 270.65
        P = 110.9063 * np.exp(-34.1632 * (h - 47) / 270.65)
    elif 51000 <= h < 71000:
        T = 413.45 - 2.8 * h  
        P = 66.93887 * (270.65 / (270.65 - 2.8 * (h - 51))) ** (34.1632 / -2.8)  
    elif 71000 <= h < 85000:
        T = 356.65 - 2.0 * h
        P = 3.956420 * (214.65 / (214.65 - 2 * (h - 71))) ** (34.1632 / -2)
    elif 85000 <= h < 91000:
        b = -3.322622E-06
        c = 9.111460E-04
        d = -0.2609971
        e = 5.944694
    elif 91000 <= h < 100000:
        b = 2.873405E-05
        c = -0.008492037
        d = 0.6541179
        e = -23.62010
    elif 100000 <= h < 110000:
        b = 0.005162063
        c = -0.8048342
        d = 55.55996
        e = -1443.338
    elif 110000 <= h < 120000:
        b = -8.854164E-05
        c = 0.03373254
        d = -4.390837
        e = 176.5294

    if 32000 < h < 85000:
        return P/(R * T)
    elif 85000 < h < 120000:
        r0 = 6356.766
        z = (-h * r0)/(h - r0)
        return np.exp(b * z ** 3 + c * z ** 2 + d * z + e)
    else:
        return 0