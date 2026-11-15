#Name: Judy Zhu
#Collaborator: Erina L.
#all angles are in radians
import numpy as np
e = 0.6615482135209798
a = 2.515867717273818 
I = 0.05948386
Omega = 4.73622346
omega = 2.34504823
M_epoch = 0.70515797
t_epoch = 2460800.5 
t_mid = 2460860.5
GM = 0.0002959151762 #GM [units] = AU^3 / day^2
epsilon = 0.40903466 #tilt of earth
n = (GM / a**(3))**0.5
R = np.array([-2.103885931031447E-01, 9.946361059061014E-01, -5.964582748524840E-05]) 
def mean_anomaly(n, M_epoch, t_mid, t_epoch):
    M_mid = M_epoch + n*(t_mid - t_epoch)
    return M_mid
def eccentric_anomaly(M_mid, e):
    M_mid = mean_anomaly(n, M_epoch, t_mid, t_epoch)
    E0 = 3
    for i in range(10000):
        E = E0 - ((E0 - e*np.sin(E0)-M_mid)/(1-e*np.cos(E0)))
        E0 = E
    return E
def position_asteroid(a, e):
    M_mid = mean_anomaly(n, M_epoch, t_mid, t_epoch)
    E = eccentric_anomaly(M_mid, e)
    pos_x = a*(np.cos(E) - e)
    pos_y = a*np.sqrt(1-e**2)*np.sin(E)
    r = np.sqrt(pos_x**2 + pos_y**2)
    return r
def true_anomaly(r):
    M_mid = mean_anomaly(n, M_epoch, t_mid, t_epoch)
    E = eccentric_anomaly(M_mid, e)
    pos_x = a*(np.cos(E) - e)
    pos_y = a*np.sqrt(1-e**2)*np.sin(E)
    r = position_asteroid(a, e)
    cos_f = pos_x / r
    sin_f = pos_y / r
    f = np.arctan2(sin_f, cos_f)
    return f
def ecliptic_r_coordinates(r, f, omega, Omega, I):
    r = position_asteroid(a, e)
    f = true_anomaly(r)
    pos_X = r*(np.cos(f+omega)*np.cos(Omega) - np.cos(I)*np.sin(f+omega)*np.sin(Omega))
    pos_Y = r*(np.cos(I)*np.cos(Omega)*np.sin(f+omega)+np.cos(f+omega)*np.sin(Omega))
    pos_Z = r*np.sin(I)*np.sin(f+omega)
    r_ecliptic = np.array([pos_X, pos_Y, pos_Z])
    return r_ecliptic
def earth_to_asteroid(R, r_ecliptic, r, f, epsilon):
    r = position_asteroid(a, e)
    f = true_anomaly(r)
    r_ecliptic = ecliptic_r_coordinates(r, f, omega, Omega, I)
    rho_X = R[0] + r_ecliptic[0] 
    rho_Y = R[1] + r_ecliptic[1]
    rho_Z = R[2] + r_ecliptic[2]
    # equatorial coordinates
    rho_x = rho_X
    rho_y = rho_Y * np.cos(epsilon) - rho_Z * np.sin(epsilon)
    rho_z = rho_Y * np.sin(epsilon) + rho_Z * np.cos(epsilon)
    rho_equatorial_vec = np.array([rho_x, rho_y, rho_z])
    return rho_equatorial_vec
def declination(r, rho_equatorial_vec, r_ecliptic, f):
    r = position_asteroid(a, e)
    f = true_anomaly(r)
    r_ecliptic = ecliptic_r_coordinates(r, f, omega, Omega, I)
    rho_equatorial_vec = earth_to_asteroid(R, r_ecliptic, r, f, epsilon)
    rho_equatorial_mag = (np.sum(rho_equatorial_vec**2))**0.5
    dec = np.rad2deg(np.arcsin(rho_equatorial_vec[2] / rho_equatorial_mag))
    return dec
def right_ascension(r, f, rho_equatorial_vec, dec, r_ecliptic):
    r = position_asteroid(a, e)
    f = true_anomaly(r)
    r_ecliptic = ecliptic_r_coordinates(r, f, omega, Omega, I)
    rho_equatorial_vec = earth_to_asteroid(R, r_ecliptic, r, f, epsilon)
    rho_equatorial_mag = declination(r, rho_equatorial_vec, r_ecliptic, f)
    cos_ra = rho_equatorial_vec[0] / (rho_equatorial_mag*np.cos(dec))
    sin_ra = rho_equatorial_vec[1] / (rho_equatorial_mag*np.cos(dec))
    ra = np.rad2deg(np.arctan2(sin_ra, cos_ra))
    return ra
def deg_to_dms(r, f, r_ecliptic, rho_equatorial_vec, dec):
    r = position_asteroid(a, e)
    f = true_anomaly(r)
    r_ecliptic = ecliptic_r_coordinates(r, f, omega, Omega, I)
    rho_equatorial_vec = earth_to_asteroid(R, r_ecliptic, r, f, epsilon)
    dec = declination(r, rho_equatorial_vec, r_ecliptic, f)
    deg = int(dec)
    minutes = int((dec - deg)*60)
    sec = ((dec - deg)*60 - minutes)*60
    return deg, minutes, sec
def deg_to_hms(r, f, rho_equatorial_vec, dec, r_ecliptic, ra):
    r = position_asteroid(a, e)
    f = true_anomaly(r)
    r_ecliptic = ecliptic_r_coordinates(r, f, omega, Omega, I)
    rho_equatorial_vec = earth_to_asteroid(R, r_ecliptic, r, f, epsilon)
    rho_equatorial_mag = declination(r, rho_equatorial_vec, r_ecliptic, f)
    ra = right_ascension(r, f, rho_equatorial_vec, dec, r_ecliptic)
    hours = int((ra) / 15)
    minutes = int(((ra / 15) - hours)*60)
    sec = (((ra / 15) - hours)*60 - minutes)*60
    return hours, minutes, sec   
def main():
    M_mid = mean_anomaly(n, M_epoch, t_mid, t_epoch)
    E = eccentric_anomaly(M_mid, e)
    r = position_asteroid(a, e)
    f = true_anomaly(r)
    r_ecliptic = ecliptic_r_coordinates(r, f, omega, Omega, I)
    rho_equatorial_vec = earth_to_asteroid(R, r_ecliptic, r, f, epsilon)
    dec = declination(r, rho_equatorial_vec, r_ecliptic, f)
    ra = right_ascension(r, f, rho_equatorial_vec, dec, r_ecliptic)
    dec_dms = deg_to_dms(r, f, r_ecliptic, rho_equatorial_vec, dec)
    ra_hms = deg_to_hms(r, f, rho_equatorial_vec, dec, r_ecliptic, ra)
    return M_mid, E, r, f, r_ecliptic, rho_equatorial_vec, dec, ra, dec_dms, ra_hms
    

    
    
    
    