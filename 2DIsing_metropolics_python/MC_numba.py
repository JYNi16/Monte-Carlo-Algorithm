# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 09:43:02 2022

@author: 26526
"""
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool 
from numba import jit 
import time
import os

sweeps = 10000
K = 1
J = 1
H = 0
Lattice = 40
relax = 100

@jit
def metropolis(spin, t):
    mag = [] 
    mag_2 = []
    for k in range(sweeps + relax):
        for i in range(Lattice):
            for j in range(Lattice):
                tmp = []
                for kk in range(8):
                    tmp.append(np.random.randint(0, Lattice))
                x_site = tmp[1]
                y_site = tmp[4]
                top = (x_site - 1) % Lattice
                bottom = (x_site + 1) % Lattice
                left = (y_site - 1) % Lattice
                right = (y_site + 1) % Lattice
                delta_E = 2.0 * spin[x_site,y_site] * (J * (spin[top,y_site] + spin[bottom,y_site] + spin[x_site,left] + spin[x_site,right]) +  H )
                if delta_E <= 0:
                    spin[x_site,y_site] = -spin[x_site,y_site]
                else:
                    if np.exp((-delta_E)/t) > np.random.random():
                        spin[x_site,y_site] = -spin[x_site,y_site]
                #spin[i,j] = 1
        
        M = np.abs(np.mean(spin))
        #print("M is:", M)
        M_2 = np.square(M)
        mag_2.append(M_2)
        mag.append(M)
    
    return mag, mag_2

        

def main_loop(T):
    
    t = T
    spin = np.ones([Lattice, Lattice], dtype = int)
    
    mag, mag_2 = metropolis(spin, t)
    
    M_ave = np.mean(mag)
    print("t is", t, "M_ave is:", M_ave)
    
    
    # print("result is:", spin)
    M_2 = np.mean(mag_2)
    C_v = (M_2 - M_ave ** 2) / T  # 求磁比热
    
    return M_2
    
def save_data(data):
    save_path = "M_T"
    if os.path.exists(save_path):
        print("save path {} exist".format(save_path))
    else:
        print("save path {} not exist".format(save_path))
        os.makedirs(save_path)
        print("now makefir the save_path")

    np.savetxt(save_path + "/M_T.txt", data)



def single_run():
    T = [round(t,2) for t in np.linspace(0.8, 4, 33)]
    M = []
    for i in range(len(T)):
        M.append(main_loop(T[i]))
        
    data = [] 
    print("M is:", M)
    for i in range(len(T)):
        data.append([T[i], M[i]])
        
    save_data(data)


def multi_run():
    start = time.time()
        
    M = []
    ing_argv = [] 
    T = [round(t,2) for t in np.linspace(0.8, 4, 33)]
    for i in range(len(T)):
        ing_argv.append(T[i])
        
    with Pool(8) as p:
        M.append(p.map(main_loop, ing_argv))
        
    print(M[0])
        
    data = [] 
    for i in range(len(T)):
        data.append([T[i], M[0][i]])
        
    save_data(data)

if __name__=="__main__":
    single_run()
