# Copyright (c) 2021 Vrund Patel, USA
#
import numpy as np
import planetary_data as pd

def coe_to_cart(a, e = 0, i = 0, RAAN = 0, argP = 0, theta = 0, u = 0, w = 0, l = 0, mu = pd.earth["mu"]):
    
    if e < 1e-11 and i != 0:
        
        e = 0
        theta = None
        argP = None
        
        r_vec_pqw = [(a * (1 - e**2)) * np.cos(u) / (1 + e * np.cos(u)), (a * (1 - e**2)) * np.sin(u) / (1 + e * np.cos(u)), 0]
        v_vec_pqw = [np.sqrt(mu / (a * (1 - e**2))) * -np.sin(u), np.sqrt(mu / (a * (1 - e**2))) * (e + np.cos(u)), 0]
        
        R1 = [[1, 0, 0],
             [0, np.cos(-i), np.sin(-i)],
             [0, -np.sin(-i), np.cos(-i)]]
        
        R3 = [[np.cos(-RAAN), np.sin(-RAAN), 0],
             [-np.sin(-RAAN), np.cos(-RAAN), 0],
             [0, 0, 1]]
        
        Q = np.dot(R3, R1)
        
        r_vec_ijk = np.dot(Q, r_vec_pqw)
        v_vec_ijk = np.dot(Q, v_vec_pqw)
        
        cart_array = np.zeros((1,6))
        cart_array[0, :3] = r_vec_ijk[:]
        cart_array[0,  3:] = v_vec_ijk[:]
        
        return cart_array
    
    elif i < 1e-11 and e != 0:
        
        i = 0
        argP = None
        RAAN = None
        
        r_vec_pqw = [(a * (1 - e**2)) * np.cos(theta) / (1 + e * np.cos(theta)), (a * (1 - e**2)) * np.sin(theta) / (1 + e * np.cos(theta)), 0]
        v_vec_pqw = [np.sqrt(mu / (a * (1 - e**2))) * -np.sin(theta), np.sqrt(mu / (a* (1 - e**2))) * (e + np.cos(theta)), 0]
        
        Q = [[np.cos(-w), np.sin(-w), 0],
            [-np.sin(-w), np.cos(-w), 0],
            [0, 0, 1]]
        
        r_vec_ijk = np.dot(Q, r_vec_pqw)
        v_vec_ijk = np.dot(Q, v_vec_pqw)
        
        return np.array([r_vec_ijk, v_vec_ijk])
    
    elif e < 1e-11 and i < 1e-11:
        
        e = 0
        i = 0
        theta = None
        argP = None
        RAAN = None
        
        r_vec_pqw = [(a * (1 - e**2)) * np.cos(l) / (1 + e*np.cos(l)), (a * (1 - e**2)) * np.sin(l) / (1 + e*np.cos(l)), 0]
        v_vec_pqw = [np.sqrt(mu / (a * (1 - e**2))) * -np.sin(l), np.sqrt(mu / (a * (1 - e**2))) * (e + np.cos(l)), 0]
        
        Q = [[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]]
        
        r_vec_ijk = np.dot(Q, r_vec_pqw)
        v_vec_ijk = np.dot(Q, v_vec_pqw)
        
        return np.array([r_vec_ijk, v_vec_ijk])
    
    else:
        
        r_vec_pqw = [(a * (1 - e**2)) * np.cos(theta) / (1 + e * np.cos(theta)), (a * (1 - e**2)) * np.sin(theta) / (1 + e * np.cos(theta)), 0]
        v_vec_pqw = [np.sqrt(mu / (a * (1 - e**2))) * -np.sin(theta), np.sqrt(mu / (a * (1 - e**2))) * (e + np.cos(theta)), 0]
        
        Q = [[np.cos(-argP) * np.cos(-RAAN) - np.cos(-i) * np.sin(-argP) * np.sin(-RAAN), np.cos(-RAAN) * np.sin(-argP) + np.cos(-argP) * np.cos(-i) *np.sin(-RAAN), np.sin(-RAAN) * np.sin(-i)],
            [-np.cos(-argP) * np.sin(-RAAN) - np.cos(-RAAN) * np.cos(-i) * np.sin(-argP), np.cos(-argP) * np.cos(-RAAN) * np.cos(-i) - np.sin(-argP) * np.sin(-RAAN), np.cos(-RAAN) * np.sin(-i)],
            [np.sin(-argP) * np.sin(-i), -np.cos(-argP) * np.sin(-i), np.cos(-i)]]
        
        r_vec_ijk = np.dot(Q, r_vec_pqw)
        v_vec_ijk = np.dot(Q, v_vec_pqw)
        
        return np.array([r_vec_ijk, v_vec_ijk])
    
    return None



def cart_to_coe(r_vec, v_vec, mu=398600):
    
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    specific_energy = (v**2 / 2) - mu / r
    
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    
    e_vec = (1 / mu) * (np.dot((v**2 - mu / r),r_vec) - np.dot((np.dot(r_vec, v_vec)), v_vec))
    e = np.linalg.norm(e_vec)
    
    a = -mu / (2 * specific_energy)
    
    i = np.rad2deg(np.arccos(h_vec[2] / h))
    
    if e < 1e-11 and i != 0:
        
        argP = None
        theta = None
        w = None
        l = None
        
        n_vec = np.cross([0, 0, 1], h_vec)
        n = np.linalg.norm(n_vec)
        RAAN = np.rad2deg(np.arccos(n_vec[0] / n))
        
        n_hat = np.dot(1 / n, n_vec)
        u = np.rad2deg(np.arccos(dot(n_hat, r_vec) / r))
        if r_vec[2] < 0:
            u = 360 - u
        
        return [a, e, i, RAAN, argP, theta, u, w, l]
    
    elif i < 1e-11 and e != 0:
        
        RAAN = None
        argP = None
        u = None
        l = None
        
        theta = np.rad2deg(np.arccos(np.dot(r_vec, e_vec) / (r * e)))
        if np.dot(r_vec, v_vec) < 0:
            theta = 360 - theta
        
        w = np.rad2deg(np.arccos(np.dot(e_vec, [1, 0, 0]) / e))
        if e_vec[1] < 0:
            w = 360 - w
        
        return [a, e, i, RAAN, argP, theta, u, w, l]
        
    elif e < 1e-11 and i < 1e-11:
        
        theta = None
        RAAN = None
        argP = None 
        u = None
        w = None
        
        l = np.rad2deg(np.arccos(np.dot(r_vec, [1, 0, 0]) / r))
        if r_vec[1] < 0:
            l = 360 - l
        
        return [a, e, i, RAAN, argP, theta, u, w, l]
    
    else:
        
        u = None
        w = None
        l = None
        
        n_vec = np.cross([0, 0, 1], h_vec)
        n = np.linalg.norm(n_vec)
        RAAN = np.rad2deg(np.arccos(n_vec[0] / n))
    
        argP = np.rad2deg(np.arccos(np.dot(n_vec, e_vec) / (n * e)))
    
        theta = np.rad2deg(np.arccos(np.dot(r_vec, e_vec) / (r * e)))
        if np.dot(r_vec, v_vec) < 0:
            theta = 360 - theta
        
        return [a, e, i, RAAN, argP, theta, u, w, l]

    return None


