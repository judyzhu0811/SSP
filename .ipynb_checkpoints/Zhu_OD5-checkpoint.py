#Name: Judy Zhu
#Collaborators: Anna K., Guthrie L., Erina L., Chris X., Arjun B., Kailin K., Isabel Y.

import numpy as np
import pandas as pd
from importlib import reload
import Zhu_odlib as od
reload(od)

#Define constants
k = 0.0172020989484 #AU^3/2 / day
cAU = 173.144643267 # speed of light in AU/(mean solar)day
eps_deg = 23.4384695831 # Earth's obliquity in degrees
eps_rad = np.deg2rad(23.4384695831)
mog_test_case = pd.read_csv('mog_test_cases.csv')

def main():
    for i in range(2):
        #get the times
        t1 = mog_test_case.loc[i, 't1_jd']
        t2 = mog_test_case.loc[i, 't2_jd']
        t3 = mog_test_case.loc[i, 't3_jd']

        print("initial t values", t1, t2, t3)
        
        t1_const = mog_test_case.loc[i, 't1_jd']
        t2_const = mog_test_case.loc[i, 't2_jd']
        t3_const = mog_test_case.loc[i, 't3_jd']
        
        t0 = mog_test_case.loc[i, 't_target_jd']
    
        #get right ascensions
        ra1 = (mog_test_case.loc[i, 'ra1']).split()
        ra1 = np.array([float(item) for item in ra1])
        ra1_rad = od.hms_to_rad(ra1)
        
        ra2 = (mog_test_case.loc[i, 'ra2']).split()
        ra2 = np.array([float(item) for item in ra2])
        ra2_rad = od.hms_to_rad(ra2)
        
        ra3 = (mog_test_case.loc[i, 'ra3']).split()
        ra3 = np.array([float(item) for item in ra3])
        ra3_rad = od.hms_to_rad(ra3)
        
        #get declinations
        dec1 = (mog_test_case.loc[i, 'dec1']).split()
        dec1 = np.array([float(item) for item in dec1])
        dec1_rad = od.dms_to_rad(dec1)
        
        dec2 = (mog_test_case.loc[i, 'dec2']).split()
        dec2 = np.array([float(item) for item in dec2])
        dec2_rad = od.dms_to_rad(dec2)
        
        dec3 = (mog_test_case.loc[i, 'dec3']).split()
        dec3 = np.array([float(item) for item in dec3])
        dec3_rad = od.dms_to_rad(dec3)

        print("RA_rad:", ra1_rad, ra2_rad, ra3_rad)
        print("DEC_rad:", dec1_rad, dec2_rad, dec3_rad)

        pos_ast_actual = np.array(mog_test_case.loc[i, 'jpl_dist2'])
        r2_vec_actual = np.array(mog_test_case.iloc[i, [20, 21, 22]])
        print("real r2_vec", r2_vec_actual)
        v2_vec_actual = np.array(mog_test_case.iloc[i, [23, 24, 25]])
        print("real v2_vec", v2_vec_actual)

        a_actual = np.array(mog_test_case.loc[i, 'jpl_a'])
        e_actual = np.array(mog_test_case.loc[i, 'jpl_e'])
        i_actual = np.array(mog_test_case.loc[i, 'jpl_i'])
        Omega_actual = np.array(mog_test_case.loc[i, 'jpl_loan'])
        omega_actual = np.array(mog_test_case.loc[i, 'jpl_aop'])
        M0_actual = np.array(mog_test_case.loc[i, 'jpl_ma'])
        
        #define rotation matrix to convert from ecliptic to equatorial
        matrix_ec2eq = np.array([[1, 0, 0],
                             [0, np.cos(eps_rad), -1*np.sin(eps_rad)],
                             [0, np.sin(eps_rad), np.cos(eps_rad)]])
        
        #define rotation matrix to convert from equatorial to ecliptic
        matrix_eq2ec = np.array([[1, 0, 0],
                             [0, np.cos(eps_rad), np.sin(eps_rad)],
                             [0, -1*np.sin(eps_rad), np.cos(eps_rad)]])
        print("matrix ec2eq:", matrix_ec2eq)
        print("matrix eq_ec:", matrix_eq2ec)
    
        #R_sun_vec in ecliptic coordinates
        R_sun_vec1 = np.array([mog_test_case.loc[i, 'sun_x_1'], mog_test_case.loc[i, 'sun_y_1'], mog_test_case.loc[i, 'sun_z_1']])
        R_sun_vec2 = np.array([mog_test_case.loc[i, 'sun_x_2'], mog_test_case.loc[i, 'sun_y_2'], mog_test_case.loc[i, 'sun_z_2']])
        R_sun_vec3 = np.array([mog_test_case.loc[i, 'sun_x_3'], mog_test_case.loc[i, 'sun_y_3'], mog_test_case.loc[i, 'sun_z_3']])

        print("R_sun vectors_ec", R_sun_vec1, R_sun_vec2, R_sun_vec3)
        #convert from ecliptic to equatorial
        R_sun_vec1_eq = matrix_ec2eq @ R_sun_vec1
        R_sun_vec2_eq = matrix_ec2eq @ R_sun_vec2
        R_sun_vec3_eq = matrix_ec2eq @ R_sun_vec3

        print("R_sun vectors_eq", R_sun_vec1_eq, R_sun_vec2_eq, R_sun_vec3_eq)
        
        #get rho_hat vectors
        rho_hat1 = np.array([np.cos(ra1_rad)*np.cos(dec1_rad), np.sin(ra1_rad)*np.cos(dec1_rad), np.sin(dec1_rad)])
        rho_hat2 = np.array([np.cos(ra2_rad)*np.cos(dec2_rad), np.sin(ra2_rad)*np.cos(dec2_rad), np.sin(dec2_rad)])
        rho_hat3 = np.array([np.cos(ra3_rad)*np.cos(dec3_rad), np.sin(ra3_rad)*np.cos(dec3_rad), np.sin(dec3_rad)])

        print("rho_hat", rho_hat1, rho_hat2, rho_hat3)
        #get tau
        tau_1 = k*(t1 - t2) #negative!!!
        tau_3 = k*(t3 - t2)
        tau_0 = k*(t3 - t1)

        print("tau values", tau_1, tau_3, tau_0)
        #constant c2 (doesn't change)
        c2 = -1
    
        #initial  for c
        c1 = tau_3/tau_0
        c3 = -tau_1/tau_0

        print("c1, c3 initial", c1, c3)
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

        print("determinants", D0, D11, D12, D13, D21, D22, D23, D31, D32, D33)
        #find rho magnitudes
        rho_1 = (c1*D11 + c2*D12 + c3*D13)/(c1*D0)
        rho_2 = (c1*D21 + c2*D22 + c3*D23)/(c2*D0)
        rho_3 = (c1*D31 + c2*D32 + c3*D33)/(c3*D0)

        print("rho_mag", rho_1, rho_2, rho_3)
        #find rho_vec
        rho1_vec = rho_1*rho_hat1
        rho2_vec = rho_2*rho_hat2
        rho3_vec = rho_3*rho_hat3

        print("rho_vec values:", rho1_vec, rho2_vec, rho3_vec)
        
        #get initial r vectors in equatorial
        r1_vec_initial = - R_sun_vec1_eq + rho1_vec
        r2_vec_initial = - R_sun_vec2_eq + rho2_vec
        r3_vec_initial = - R_sun_vec3_eq + rho3_vec
        
        r2_mag = np.sqrt(np.sum(r2_vec_initial**2))
        r2i_mag = r2_mag
        
        print("R2_initial mag:", r2i_mag)
        print("r_vec initial in equatorial", r1_vec_initial, r2_vec_initial, r3_vec_initial)
        #get initial v vectors
        v12_vec = (r2_vec_initial - r1_vec_initial)/(-1*tau_1)
        v23_vec = (r3_vec_initial - r2_vec_initial)/(tau_3)
        v2_vec = v12_vec*(tau_3/tau_0) - v23_vec*(tau_1/tau_0) #inital v2

        print("initial velocity", v12_vec, v23_vec, v2_vec)
        
        #get coefficients
        u = 1 / (r2_mag**3)
        z = np.dot(r2_vec_initial, v2_vec) / (r2_mag**2)
        q = np.dot(v2_vec, v2_vec) / (r2_mag**2) - u

        print("u, z, q", u, z, q)
        #get initial f and g values
        f1 = 1 - 0.5*u*tau_1**2 + 0.5*u*z*tau_1**3 + (1/24)*(3*u*q-15*u*z**2+u**2)*tau_1**4
        g1 = tau_1 - (1/6)*u*tau_1**3 + (1/4)*u*z*tau_1**4
        f3 = 1 - 0.5*u*tau_3**2 + 0.5*u*z*tau_3**3 + (1/24)*(3*u*q-15*u*z**2+u**2)*tau_3**4
        g3 = tau_3 - (1/6)*u*tau_3**3 + (1/4)*u*z*tau_3**4

        print("initial fg values", f1, g1, f3, g3)
        #calculate new initial c and d values
        c1 = g3 / (f1*g3 - g1*f3)
        c3 = -1*g1 / (f1*g3 - g1*f3)
        d1 = -1*f3 / (f1*g3 - g1*f3)
        d3 = f1 / (f1*g3 - g1*f3)

        print("initial cd values", c1, c3, d1, d3)
        #calculate new r2 and v2
        r2_vec = c1*r1_vec_initial + c3*r3_vec_initial
        v2_vec = d1*r1_vec_initial + d3*r3_vec_initial
        r2f_mag = np.sqrt(np.sum(r2_vec**2))

        print("first iteration r2_vec and v2_vec", r2_vec, v2_vec)
        print("first iteration r2 mag:", r2f_mag)
        
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
            r2f_mag = r2_mag
            
            #v2_vec in equatorial, AU/days
            v2_vec = d1*r1_vec + d3*r3_vec
            
            #convert from equatorial to ecliptic coordinates
            r1_vec_ec = matrix_eq2ec @ r1_vec
            r2_vec_ec = matrix_eq2ec @ r2_vec
            r3_vec_ec = matrix_eq2ec @ r3_vec

            #get magnitude of r2 in ecliptic
            r2_mag_ec = np.linalg.norm(r2_vec_ec)
            
            #v2_vec in ecliptic and Gaussian days
            v2_vec_ec = k*(d1*r1_vec_ec + d3*r3_vec_ec)
              
            count += 1
            if count ==2000:
                print("didnt converge :((")
                break
                
        # asteroid position in ecliptic
        print("equatorial", r2_vec, v2_vec)
        pos_ast_vec = matrix_eq2ec @ rho2_vec
        pos_ast_mag = np.linalg.norm(pos_ast_vec)
        pos_ast_error = 100* abs(pos_ast_mag - np.linalg.norm(pos_ast_actual)) / np.linalg.norm(pos_ast_actual)
        
# get r2_vec_ec, v2_vec_ec, fg values, and % errors
        pos_error = 100* abs(r2_mag_ec - np.linalg.norm(r2_vec_actual)) / np.linalg.norm(r2_vec_actual)
        vel_error = 100* abs(np.linalg.norm(v2_vec_ec) - np.linalg.norm(v2_vec_actual)) / np.linalg.norm(v2_vec_actual)
        print("c1, c3 after", c1, c3)
        print("fg values", f1, g1, f3, g3)
        print("tau values", tau_1, tau_3, tau_0)
        print("number of iterations", count)
        print("rho_values after:", rho_1, rho_2, rho_3)
        print("r2i, r2f", r2i_mag, r2f_mag)
        print("asteroid distance error:", pos_ast_error, "%")
        print("r_vec_ec", r1_vec_ec, r2_vec_ec, r3_vec_ec)
        print("v2_vec_ec", v2_vec_ec)
        print("f1, g1, f3, g3", f1, g1, f3, g3)
        print("% error r2_vec", pos_error, "%")
        print("% error v2_vec", vel_error, "%")

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
        
        M_epoch = od.mean_anomaly_epoch(r2_mag_ec, r2_vec_ec, v2_vec_ec/k, a, e, t0, t2_const, Omega, I, f)
        M_epoch_deg = np.rad2deg(M_epoch)
        print(h_vec, a, e, I, Omega, omega, M_epoch_deg)
        
        # get orbital elements % errors

        a_error = 100*abs(a - a_actual)/a_actual
        e_error = 100*abs(e - e_actual)/e_actual
        i_error = 100*abs(I_deg - i_actual)/i_actual
        Omega_error = 100*abs(Omega_deg - Omega_actual)/Omega_actual
        omega_error = 100*abs(omega_deg - omega_actual)/omega_actual
        M_epoch_error = 100*abs(M_epoch_deg - M0_actual)/M0_actual

        print("Orbital Elements error %:", a_error, e_error, i_error, Omega_error, omega_error, M_epoch_error)
    
            


            
            
            
