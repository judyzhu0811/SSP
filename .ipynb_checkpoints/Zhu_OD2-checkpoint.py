#Name: Judy Zhu
#Date: 7/18/26
#Collaborators: Chris, Phong
import numpy as np
pos_vec = np.array([0.244, 2.17, -0.445]) #AU/Gaussian day
pos_mag = np.sum(pos_vec**2)**0.5 #Au/Gaussian day
vel_vec = np.array([-0.731, -0.0041, 0.0502]) #AU/Gaussian day
h_cross = np.cross(pos_vec, vel_vec)
h_mag = (np.sum(np.cross(pos_vec, vel_vec)**2))**0.5
GM = 0.000295912208231 #GM [units] = AU^3 / day^2
t_epoch = 2451545.0 #JHD
t_mid = 2456842.5 #JHD

#all angles are in radians
def semimajor_axis(pos_mag, vel_vec):
    
    a = ( 2/pos_mag - np.dot(vel_vec, vel_vec))**(-1)
    return a
    
def eccentricity(pos_mag, vel_vec, a, h_mag):
    
    a = semimajor_axis(pos_mag, vel_vec)
    e = (1 - h_mag**2 / a) ** 0.5
    return e
    
def inclination(h_cross, h_mag):
    I = np.arccos(h_cross[2] / h_mag) # angle in radians
    return I
    
def longitude_ascending_node(h_cross):
    
    Omega = np.arctan2(h_cross[0], -1 * h_cross[1])
    return Omega
    
def longitude_perihelion(pos_vec, pos_mag, vel_vec, Omega, a, e, I):
    
    #definitions
    a = semimajor_axis(pos_mag, vel_vec)
    e = eccentricity(pos_mag, vel_vec, a, h_mag)
    I = inclination(h_cross, h_mag)
    Omega = longitude_ascending_node(h_cross)
    
    #get f + omega
    sin_total = pos_vec[2] / (pos_mag * np.sin(I))
    cos_total = ((pos_vec[0] / pos_mag) + np.cos(I) * (pos_vec[2] / (pos_mag * np.sin(I)) * np.sin(Omega))) * (1 / np.cos(Omega))
    f_Omega = np.arctan2(sin_total, cos_total)
    
    #get f
    sin_f = np.dot(pos_vec, vel_vec) * ((a * (1 - e**2)) **0.5) / (e*pos_mag)
    cos_f = (a*(1-e**2)/pos_mag - 1)/e
    f = np.arctan2(sin_f, cos_f)
    
    #get omega
    omega = f_Omega - f
    return omega 
    
def mean_anomaly_epoch(pos_mag, a, e, t_epoch, t_mid, Omega, I):
    
    f = longitude_perihelion(pos_vec, pos_mag, vel_vec, Omega, a, e, I)
    if f>0:
        E = np.arccos((1 - pos_mag / a) / e) 
    else:
        E = -1 * np.arccos((1 - pos_mag / a) / e)
    n = (GM / a**(3))**0.5
    M_mid = E - e*(np.sin(E))
    M_epoch = (M_mid + n*(t_epoch - t_mid))%(np.pi*2)
    return M_epoch
    
def main():  
    a = semimajor_axis(pos_mag, vel_vec)
    e = eccentricity(pos_mag, vel_vec, a, h_mag)
    I = inclination(h_cross, h_mag)
    Omega = longitude_ascending_node(h_cross)
    omega = longitude_perihelion(pos_vec, pos_mag, vel_vec, Omega, a, e, I) + 2*np.pi
    M_epoch = mean_anomaly_epoch(pos_mag, a, e, t_epoch, t_mid, Omega, I)
    return a, e, I, Omega, omega, M_epoch
    
