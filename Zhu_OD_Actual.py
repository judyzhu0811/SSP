#Name: Judy Zhu
import numpy as np
import pandas as pd
from importlib import reload
import Zhu_odlib as od
reload(od)

#define rotation matrix to convert from ecliptic to equatorial
image_data = pd.read_csv('IMAGE_DATA.csv')

#Define constants
k = 0.0172020989484 #AU^3/2 / day
cAU = 173.144643267 # speed of light in AU/(mean solar)day
eps_deg = 23.4384695831 # Earth's obliquity in degrees
eps_rad = np.deg2rad(23.4384695831)

#rotation matrix from ecliptic to equatorial
matrix_ec2eq = np.array([[1, 0, 0],
                    [0, np.cos(eps_rad), -1*np.sin(eps_rad)],
                    [0, np.sin(eps_rad), np.cos(eps_rad)]])
        
#define rotation matrix to convert from equatorial to ecliptic
matrix_eq2ec = np.array([[1, 0, 0],
                    [0, np.cos(eps_rad), np.sin(eps_rad)],
                    [0, -1*np.sin(eps_rad), np.cos(eps_rad)]])
#get the times (this is halfway)
t1 = 2461227.6979167 
t2 = 2461231.7500000 
t3 = 2461236.7083333
t4 = 2461240.6666667

t1_const = 2461227.6979167 
t2_const = 2461231.7500000 
t3_const = 2461236.7083333
t0 = 2458473.5 
t4_const = 2461240.6666667

#times from pics
# t1 = 2461227.7055556
# t2 = 2461231.7388889
# t3 = 2461236.7034722

# t1_const = 2461227.7055556
# t2_const = 2461231.7388889
# t3_const = 2461236.7034722
# t0 = 2461231.7388889

#Plate solving angles from our observations
# ra1_rad = np.deg2rad(246.864008)
# dec1_rad = np.deg2rad(8.96303333)
# ra2_rad = np.deg2rad(245.633742)
# dec2_rad = np.deg2rad(10.8217722)
# ra3_rad = np.deg2rad(244.40145)
# dec3_rad = np.deg2rad(12.8442111)
# ra4_rad = np.deg2rad(243.684262)
# dec4_rad = np.deg2rad(14.2092056)

#Use JPL Angles
ra1_rad = np.deg2rad(246.866541667)
dec1_rad = np.deg2rad(8.95958333333)
ra2_rad = np.deg2rad(245.6305)
dec2_rad = np.deg2rad(10.8265277778)
ra3_rad = np.deg2rad(244.412666667)
dec3_rad = np.deg2rad(12.8259166667)

#get R_sun_vec in ecliptic coordinates values in AU/day
r2_vec_actual = np.array([6.597708943023856E-02, -1.414983210141008E+00, 3.160210310129223E-01])
v2_vec_actual = np.array([1.288134902030575E-02, 6.902634010416600E-03, 4.594750177786959E-03])
pos_ast_actual = np.array([-2.423378396237330E-01, -4.462283685152122E-01, 3.159296182543370E-01])

#get orbital elements in degrees
a_actual = 1.709375875200247    
e_actual = 0.3648329026171276   
i_actual = 26.86986395245218
Omega_actual = 246.5733908423288   
omega_actual = 115.4656355507939   
M0_actual = 3.134881025916907E+02

#Units are in AU and days (from observations)
R_sun_vec1 = np.array([-2.423315750408254E-01, 9.873394152864740E-01, -5.667175240168493E-05])
R_sun_vec2 = np.array([-3.083053310421857E-01, 9.687376207018604E-01, -5.367674882323271E-05])
R_sun_vec3 = np.array([-3.870854813412044E-01, 9.398635835173186E-01, -5.292927887796023E-05])
R_sun_vec4 = np.array([-4.486973102323266E-01, 9.117911989410320E-01, -5.346925197728322E-05])

#get R_sun_vec in ecliptic coordinates values in AU/day
r2_vec_actual = np.array([6.597708943023856E-02, -1.414983210141008E+00, 3.160210310129223E-01])
v2_vec_actual = np.array([1.288134902030575E-02, 6.902634010416600E-03, 4.594750177786959E-03])
pos_ast_actual = np.array([-2.423378396237330E-01, -4.462283685152122E-01, 3.159296182543370E-01])

#get orbital elements in degrees
a_actual = 1.709375875200247    
e_actual = 0.3648329026171276   
i_actual = 26.86986395245218
Omega_actual = 246.5733908423288   
omega_actual = 115.4656355507939   
M0_actual = 3.134881025916907E+02

#Units are in AU and days (from observations)
R_sun_vec1 = np.array([-2.423315750408254E-01, 9.873394152864740E-01, -5.667175240168493E-05])
R_sun_vec2 = np.array([-3.083053310421857E-01, 9.687376207018604E-01, -5.367674882323271E-05])
R_sun_vec3 = np.array([-3.870854813412044E-01, 9.398635835173186E-01, -5.292927887796023E-05])
R_sun_vec4 = np.array([-4.486973102323266E-01, 9.117911989410320E-01, -5.346925197728322E-05])

def get_od_elements(t1, t2, t3, t1_const, t2_const, t3_const, ra1_rad, dec1_rad, ra2_rad, dec2_rad, ra3_rad, dec3_rad, r2_vec_actual, v2_vec_actual, a_actual, e_actual, i_actual, Omega_actual, omega_actual, M0_actual, R_sun_vec1, R_sun_vec2, R_sun_vec3):
    

    #convert from ecliptic to equatorial
    R_sun_vec1_eq = matrix_ec2eq @ R_sun_vec1
    R_sun_vec2_eq = matrix_ec2eq @ R_sun_vec2
    R_sun_vec3_eq = matrix_ec2eq @ R_sun_vec3
    
    
    #get rho_hat vectors
    rho_hat1 = np.array([np.cos(ra1_rad)*np.cos(dec1_rad), np.sin(ra1_rad)*np.cos(dec1_rad), np.sin(dec1_rad)])
    rho_hat2 = np.array([np.cos(ra2_rad)*np.cos(dec2_rad), np.sin(ra2_rad)*np.cos(dec2_rad), np.sin(dec2_rad)])
    rho_hat3 = np.array([np.cos(ra3_rad)*np.cos(dec3_rad), np.sin(ra3_rad)*np.cos(dec3_rad), np.sin(dec3_rad)])
    
    #get tau
    tau_1 = k*(t1 - t2) #negative!!!
    tau_3 = k*(t3 - t2)
    tau_0 = k*(t3 - t1)

    #initial  for c
    c2 = -1
    c1 = tau_3/tau_0
    c3 = -tau_1/tau_0

    #determinants
    D0 = np.dot(rho_hat1, np.cross(rho_hat2, rho_hat3))
    D11 = np.dot(np.cross(R_sun_vec1_eq, rho_hat2), rho_hat3)
    D12 = np.dot(np.cross(R_sun_vec2_eq, rho_hat2), rho_hat3)
    D13 = np.dot(np.cross(R_sun_vec3_eq, rho_hat2), rho_hat3)
    D21 = np.dot(np.cross(rho_hat1, R_sun_vec1_eq), rho_hat3)
    D22 = np.dot(np.cross(rho_hat1, R_sun_vec2_eq), rho_hat3)
    D23 = np.dot(np.cross(rho_hat1, R_sun_vec3_eq), rho_hat3)
    D31 = np.dot(rho_hat1, np.cross(rho_hat2, R_sun_vec1_eq))
    D32 = np.dot(rho_hat1, np.cross(rho_hat2, R_sun_vec2_eq))
    D33 = np.dot(rho_hat1, np.cross(rho_hat2, R_sun_vec3_eq))

    #find rho magnitudes
    rho_1 = (c1*D11 + c2*D12 + c3*D13)/(c1*D0)
    rho_2 = (c1*D21 + c2*D22 + c3*D23)/(c2*D0)
    rho_3 = (c1*D31 + c2*D32 + c3*D33)/(c3*D0)

    #find rho_vec
    rho1_vec = rho_1*rho_hat1
    rho2_vec = rho_2*rho_hat2
    rho3_vec = rho_3*rho_hat3

    
    #get initial r vectors in equatorial
    r1_vec_initial = - R_sun_vec1_eq + rho1_vec
    r2_vec_initial = - R_sun_vec2_eq + rho2_vec
    r3_vec_initial = - R_sun_vec3_eq + rho3_vec
    
    r2_mag = np.sqrt(np.sum(r2_vec_initial**2))
    r2i_mag = r2_mag

    #get initial v vectors
    v12_vec = (r2_vec_initial - r1_vec_initial)/(-1*tau_1)
    v23_vec = (r3_vec_initial - r2_vec_initial)/(tau_3)
    v2_vec = v12_vec*(tau_3/tau_0) - v23_vec*(tau_1/tau_0) #inital v2

    
    #get coefficients
    u = 1 / (r2_mag**3)
    z = np.dot(r2_vec_initial, v2_vec) / (r2_mag**2)
    q = np.dot(v2_vec, v2_vec) / (r2_mag**2) - u

    #get initial f and g values
    f1 = 1 - 0.5*u*tau_1**2 + 0.5*u*z*tau_1**3 + (1/24)*(3*u*q-15*u*z**2+u**2)*tau_1**4
    g1 = tau_1 - (1/6)*u*tau_1**3 + (1/4)*u*z*tau_1**4
    f3 = 1 - 0.5*u*tau_3**2 + 0.5*u*z*tau_3**3 + (1/24)*(3*u*q-15*u*z**2+u**2)*tau_3**4
    g3 = tau_3 - (1/6)*u*tau_3**3 + (1/4)*u*z*tau_3**4

    #calculate new initial c and d values
    c1 = g3 / (f1*g3 - g1*f3)
    c3 = -1*g1 / (f1*g3 - g1*f3)
    d1 = -1*f3 / (f1*g3 - g1*f3)
    d3 = f1 / (f1*g3 - g1*f3)

    #calculate new r2 and v2
    r2_vec = c1*r1_vec_initial + c3*r3_vec_initial
    v2_vec = d1*r1_vec_initial + d3*r3_vec_initial
    r2f_mag = np.sqrt(np.sum(r2_vec**2))

    ##START ITERATIONS##
    count = 0
    while ((abs(r2f_mag-r2i_mag)/r2f_mag) > 10**(-14)):
      
        #correct for light
        t1 = t1_const - rho_1/cAU
        t2 = t2_const - rho_2/cAU
        t3 = t3_const - rho_3/cAU

        r2i_mag = r2f_mag
        
        #recalculate tau values
        tau_1 = k*(t1 - t2) #negative!!!
        tau_3 = k*(t3 - t2)
        tau_0 = k*(t3 - t1)

        #get f, g values from the new time
        f1 = od.get_fourth_f(r2_vec, v2_vec, tau_1)
        g1 = od.get_fourth_g(r2_vec, v2_vec, tau_1)
        f3 = od.get_fourth_f(r2_vec, v2_vec, tau_3)
        g3 = od.get_fourth_g(r2_vec, v2_vec, tau_3)
        
        #recalculate c and d values
        c1 = g3 / (f1*g3 - g1*f3)
        c3 = -1*g1 / (f1*g3 - g1*f3)
        d1 = -1*f3 / (f1*g3 - g1*f3)
        d3 = f1 / (f1*g3 - g1*f3)
    
        #recalculate rho values
        rho_1 = (c1*D11 + c2*D12 + c3*D13)/(c1*D0)
        rho_2 = (c1*D21 + c2*D22 + c3*D23)/(c2*D0)
        rho_3 = (c1*D31 + c2*D32 + c3*D33)/(c3*D0)
    
        #find rho_vec
        rho1_vec = rho_1*rho_hat1
        rho2_vec = rho_2*rho_hat2
        rho3_vec = rho_3*rho_hat3

        #get initial r vectors
        r1_vec = - R_sun_vec1_eq + rho1_vec
        r2_vec = - R_sun_vec2_eq + rho2_vec
        r3_vec = - R_sun_vec3_eq + rho3_vec
        r2_mag = np.sqrt(np.sum(r2_vec**2))

        #while loop
        r2f_mag = r2_mag
        
        #v2_vec in equatorial, AU/days
        v2_vec = d1*r1_vec + d3*r3_vec
        
        #convert from equatorial to ecliptic coordinates
        r1_vec_ec = matrix_eq2ec @ r1_vec
        r2_vec_ec = matrix_eq2ec @ r2_vec
        r3_vec_ec = matrix_eq2ec @ r3_vec
          
        count += 1
        if count ==2000:
            print("didnt converge :((")
            break


    #get magnitude of r2 in ecliptic
    r2_mag_ec = np.linalg.norm(r2_vec_ec)
        
    #v2_vec in ecliptic and Gaussian days
    v2_vec_ec = k*(d1*r1_vec_ec + d3*r3_vec_ec)
    
    # asteroid position in ecliptic
    pos_ast_vec = matrix_eq2ec @ rho2_vec
    pos_ast_mag = np.linalg.norm(pos_ast_vec)
    pos_ast_error = 100* abs(pos_ast_mag - np.linalg.norm(pos_ast_actual)) / np.linalg.norm(pos_ast_actual)
    
# get r2_vec_ec, v2_vec_ec, fg values, and % errors
    pos_error = 100* abs(r2_mag_ec - np.linalg.norm(r2_vec_actual)) / np.linalg.norm(r2_vec_actual)
    vel_error = 100* abs(np.linalg.norm(v2_vec_ec) - np.linalg.norm(v2_vec_actual)) / np.linalg.norm(v2_vec_actual)

    # get orbital elements
    h_vec = od.specific_angular_momentum(r2_vec_ec, v2_vec_ec/k)[0]
    h_mag = od.specific_angular_momentum(r2_vec_ec, v2_vec_ec/k)[1]
    a = od.semimajor_axis(r2_mag_ec, v2_vec_ec/k)
    e = od.eccentricity(a, h_mag, r2_vec_ec, v2_vec_ec/k)
    I = od.inclination(h_vec, r2_vec_ec, v2_vec_ec/k, h_mag)
    I_deg = np.rad2deg(od.inclination(h_vec, r2_vec_ec, v2_vec_ec/k, h_mag))

    Omega = od.longitude_ascending_node(h_vec, r2_vec_ec, v2_vec_ec/k, I, a, h_mag)
    Omega_deg = (np.rad2deg(od.longitude_ascending_node(h_vec, r2_vec_ec, v2_vec_ec/k, I, a, h_mag)))
    
    f = np.arctan2(np.dot(r2_vec_ec, v2_vec_ec/k) * ((a * (1 - e**2)) **0.5) / (e*r2_mag_ec), (a*(1-e**2)/r2_mag_ec - 1)/e)
    omega = (od.longitude_perihelion(r2_vec_ec, r2_mag_ec, v2_vec_ec/k, h_vec, h_mag, I, Omega, a, e))
    omega_deg = (np.rad2deg(od.longitude_perihelion(r2_vec_ec, r2_mag_ec, v2_vec_ec/k, h_vec, h_mag, I, Omega, a, e)))
    
    M_epoch = od.mean_anomaly_epoch(r2_mag_ec, r2_vec_ec, v2_vec_ec/k, a, e, t0, t2_const, Omega, I, f) + np.pi
    M_epoch_deg = np.rad2deg(M_epoch)
    
    
    
    # get orbital elements % errors
    a_error = 100*abs(a - a_actual)/a_actual
    e_error = 100*abs(e - e_actual)/e_actual
    i_error = 100*abs(I_deg - i_actual)/i_actual
    Omega_error = 100*abs(Omega_deg - Omega_actual)/Omega_actual
    omega_error = 100*abs(omega_deg - omega_actual)/omega_actual
    M_epoch_error = 100*abs(M_epoch_deg - M0_actual)/M0_actual

    
    print("number of iterations", count)
    print("asteroid distance error:", pos_ast_error, "%")
    print("r_vec_ec", r1_vec_ec, r2_vec_ec, r3_vec_ec)
    print("r2_vec_actual, v2_vec_actual", r2_vec_actual, v2_vec_actual)
    print("v2_vec_ec", v2_vec_ec)
    print("f1, g1, f3, g3", f1, g1, f3, g3)
    print("% error r2_vec", pos_error, "%")
    print("% error v2_vec", vel_error, "%")


    print("Orbital Elements:", a, e, I_deg, Omega_deg, omega_deg, M_epoch_deg)
    print("Expected Orbital Elements:", a_actual, e_actual, i_actual, Omega_actual, omega_actual, M0_actual)
    print("Orbital Elements error %:", a_error, e_error, i_error, Omega_error, omega_error, M_epoch_error)
    return np.array([a, e, I, Omega_deg, omega_deg, M_epoch_deg])

            


            
            
            
