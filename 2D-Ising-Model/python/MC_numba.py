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

sweeps = 2000
K = 1
J = 1
H = 0
Lattice = 40
relax = 100

@jit
def metrics_loop(spin, t):
    mag = [] 
    mag_2 = []
    for k in range(sweeps + relax):
        for i in range(Lattice):
            for j in range(Lattice):
                top = (i - 1) % Lattice
                bottom = (i + 1) % Lattice
                left = (j - 1) % Lattice
                right = (j + 1) % Lattice
                delta_E = 2.0 * spin[i,j] * (J * (spin[top,j] + spin[bottom,j] + spin[i,left] + spin[i,right]) +  H )
                if delta_E <= 0:
                    spin[i,j] = -spin[i,j]
                else:
                    if np.exp((-delta_E)/t) > np.random.random():
                        spin[i,j] = -spin[i,j]
                #spin[i,j] = 1
        
        M = np.abs(np.mean(spin))
        #print("M is:", M)
        M_2 = np.square(M)
        mag_2.append(M_2)
        mag.append(M)
    
    return mag, mag_2
       
       
        

def main_loop(T):
    
    t = T/10 
    spin = np.ones([Lattice, Lattice], dtype = int)
    
    mag, mag_2 = metrics_loop(spin, t)
    
            
    result = np.sum(mag[relax:]) / sweeps
    print("result is:", result)
    test = [] 

        #print("result is:", spin)
    result_2 = np.sum(mag_2[relax:]) / sweeps 
    result_C = (result_2 - result ** 2) / T # 求磁比热
    test.append(t)
    test.append(result)
    test.append(result_C)
        
    np.save("M_T/{:d}.npy".format(T), test)

    
def single_run():
    start = time.time()
    for T in range(3, 30, 1):
        main_loop(T)
        
    end = time.time() 
    print("run time is:", (end - start))

def multi_run():
    start = time.time()
        
    ing_argv = [] 
    for T in range(8, 40, 1):
        ing_argv.append(T)
    with Pool(16) as p:
        p.map(main_loop, ing_argv)
        
    end = time.time() 
        
    print("run time is:", (end - start))

if __name__=="__main__":
    single_run()
