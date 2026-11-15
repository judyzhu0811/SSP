import numpy as np
import pandas as pd
from importlib import reload
#Collaborator: Chris X., Erina L.
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

def main():
    f_g_test_case = pd.read_csv('f_and_g_test_cases.csv')
    for i in range(6):
        r2_vec = f_g_test_case.iloc[i, [0, 1, 2]]
        v2_vec = f_g_test_case.iloc[i, [3, 4, 5]]
        tau = f_g_test_case.iloc[i, [6]]
        f_third = np.array(get_third_f(r2_vec, v2_vec, tau))
        f_fourth = np.array(get_fourth_f(r2_vec, v2_vec, tau))
        g_third = np.array(get_third_g(r2_vec, v2_vec, tau))
        g_fourth = np.array(get_fourth_g(r2_vec, v2_vec, tau))
        print("f_third differences:", 0.995 - 0.995, -0.99769896 + 0.997699, -0.98244939 + 0.982449)
        print("g_third differences:", 0.09983333 - 0.099833, 0.49960229 - 0.499602, -0.32427006 - (-0.32427))
        print("f_fourth differences:", -0.99499323 + 0.995, 0.99769475 - 0.997699, -0.98232462 + 0.982449)
        print("g_fourth differences:", 0.09983333 - 0.099833, 0.4996236 - 0.499602, -0.32400422 - (-0.32427))


    
    