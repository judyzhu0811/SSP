import numpy as np
import pandas as pd
from importlib import reload

#-------------OD1-----------------------------------------------#
def specific_angular_momentum(pos_vec, vel_vec):
    h_vec = np.cross(pos_vec, vel_vec) #angular momentum 
    h_mag = np.linalg.norm(h_vec)
    return h_vec, h_mag

#-------------OD2-------------------------------------------------#
def semimajor_axis(pos_mag, vel_vec):
    a = ( 2/pos_mag - np.dot(vel_vec, vel_vec))**(-1)
    return a
    
def eccentricity(a, h_mag, pos_vec, vel_vec):
    e = np.sqrt(1 - h_mag**2/a)
    return e
    
def inclination(h_vec, pos_vec, vel_vec, h_mag):
    I = np.arccos(h_vec[2] / h_mag) # angle in radians
    return I
    
def longitude_ascending_node(h_vec, pos_vec, vel_vec, I, a, h_mag):
    Omega = np.arctan2(h_vec[0], -1 * h_vec[1]) % (2*np.pi)
    print(h_vec)
    return Omega 
    
def longitude_perihelion(pos_vec, pos_mag, vel_vec, h_vec, h_mag, I, Omega, a, e):
    #get f + omega
    sin_total = pos_vec[2] / (pos_mag * np.sin(I))
    cos_total = ((pos_vec[0] / pos_mag) + np.cos(I) * sin_total * np.sin(Omega)) * (1 / np.cos(Omega))
    f_Omega = np.arctan2(sin_total, cos_total) % (2*np.pi)
    print("inclination", I)
    print(sin_total, cos_total)
    print(np.rad2deg(f_Omega))
    #get f
    sin_f = np.dot(pos_vec, vel_vec) * ((a * (1 - e**2)) **0.5) / (e*pos_mag)
    cos_f = (a*(1-e**2)/pos_mag - 1)/e
    f = np.arctan2(sin_f, cos_f) % (2*np.pi)

    print(np.rad2deg(f))
    #get omega
    omega = f_Omega - f
    return omega % (2*np.pi)
    
def mean_anomaly_epoch(pos_mag, pos_vec, vel_vec, a, e, t_epoch, t_mid, Omega, I, f):
#t in Julian days
#
    sin_f = np.dot(pos_vec, vel_vec) * ((a * (1 - e**2)) **0.5) / (e*pos_mag)
    cos_f = (a*(1-e**2)/pos_mag - 1)/e
    f = np.arctan2(sin_f, cos_f)
    if f>0:
        E = np.arccos((1 - pos_mag / a) / e) 
    else:
        E = -1 * np.arccos((1 - pos_mag / a) / e)
    # n = (GM / a**(3))**0.5
    k = 0.0172020989484
    M_mid = E - e*(np.sin(E))
    M_epoch = (M_mid + k*(t_epoch - t_mid) / a**(1.5))
    return M_epoch % (2*np.pi)
    
def main():  
    a = semimajor_axis(pos_mag, vel_vec)
    e = eccentricity(pos_mag, vel_vec, a, h_mag)
    I = inclination(h_cross, h_mag)
    Omega = longitude_ascending_node(h_cross)
    omega = longitude_perihelion(pos_vec, pos_mag, vel_vec, Omega, a, e, I)
    M_epoch = mean_anomaly_epoch(pos_mag, a, e, t_epoch, t_mid, Omega, I)
    return a, e, I, Omega, omega, M_epoch
#--------------OD3-------------------------------------------------------------------
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
#-------------------OD4------------------------------------------------------
def get_third_f(r2_vec, v2_vec, tau):
    r2_mag = (np.sum(r2_vec**2))**0.5
    v2_mag = (np.sum(v2_vec**2))**0.5
    u = 1 / (r2_mag)**3
    z = np.dot(r2_vec, v2_vec) / (r2_mag)**2
    q = np.dot(v2_vec, v2_vec)/(r2_mag)**2 - u
    f_third = 1 - 0.5*u*tau**2+0.5*u*z*tau**3
    return f_third
def get_third_g(r2_vec, v2_vec, tau):
    r2_mag = (np.sum(r2_vec**2))**0.5
    v2_mag = (np.sum(v2_vec**2))**0.5
    u = 1 / (r2_mag)**3
    z = np.dot(r2_vec, v2_vec) / (r2_mag)**2
    q = np.dot(v2_vec, v2_vec)/(r2_mag)**2 - u
    g_third = tau - (1/6)*u*tau**3
    return g_third
def get_fourth_f(r2_vec, v2_vec, tau):
    r2_mag = (np.sum(r2_vec**2))**0.5
    v2_mag = (np.sum(v2_vec**2))**0.5
    u = 1 / (r2_mag)**3
    z = np.dot(r2_vec, v2_vec) / (r2_mag)**2
    q = np.dot(v2_vec, v2_vec)/(r2_mag)**2 - u
    f_fourth = 1 - 0.5*u*tau**2+0.5*u*z*tau**3 +(1/24)*(3*u*q - 15*u*z**2 + u**2)*tau**4
    return f_fourth
def get_fourth_g(r2_vec, v2_vec, tau):
    r2_mag = (np.sum(r2_vec**2))**0.5
    v2_mag = (np.sum(v2_vec**2))**0.5
    u = 1 / (r2_mag)**3
    z = np.dot(r2_vec, v2_vec) / (r2_mag)**2
    q = np.dot(v2_vec, v2_vec)/(r2_mag)**2 - u
    g_fourth = tau - (1/6)*u*tau**3+(1/4)*u*z*tau**4
    return g_fourth
#----------------OD5-------------------------------------
def dms_to_rad(dec):
    if dec[0] > 0:
        dec_rad = np.deg2rad(abs(dec[0]) + abs(dec[1]/60) + abs(dec[2]/3600))
    else:
        dec_rad = -np.deg2rad(abs(dec[0]) + abs(dec[1]/60) + abs(dec[2]/3600))
    return dec_rad %(2*np.pi)
def hms_to_rad(ra):
    ra_rad = np.deg2rad(15*abs(ra[0]) + abs(ra[1]/4) + abs(ra[2]/240))
    return ra_rad
def hms_to_deg(ra):
    ra_deg = 15*ra
    return ra_deg
    