## Author jinyang ni 
## Monte-Carlo simulation for 2D-Ising model 

import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool 
from numba import jit 

import time

class MonteCarlo():
    
    def __init__(self, relax, sweeps, K, J, H, Lattice):
        self.sweeps = sweeps
        self.relax = relax 
        self.K = K 
        self.J = J 
        self.H = H 
        self.Lattice = Lattice
    
    def cal_delta_E(self, spin, i, j):
        
        top = (i - 1) % self.Lattice
        bottom = (i + 1) % self.Lattice
        left = (j - 1) % self.Lattice
        right = (j + 1) % self.Lattice
        E_loc_former = -1.0 * spin[i,j] * (self.J * (spin[top,j] + spin[bottom,j] + spin[i,left] + spin[i,right]) +  self.H )
        E_loc_next = spin[i,j] * (self.J * (spin[top,j] + spin[bottom,j] + spin[i,left] + spin[i,right]) +  self.H )
        delta_E = E_loc_next - E_loc_former
        
        return delta_E 
    
    #@jit
    def main_loop(self, T):
        
        S = self.Lattice ** 2
        #f_M = open('M.dat','w')
        #f_C = open('C.dat','w')
            
        #for T in np.arange(0.1, 2.0, 0.05):
        t = T/10 
        spin = np.ones([self.Lattice, self.Lattice], dtype = int)
        C = 0
        C_heat = []
        mag = []
        mag_2 = []
        for sweep in range(self.sweeps + self.relax):
            for i in range(self.Lattice):
                for j in range(self.Lattice):
                    tmp = []
                    for kk in range(8):
                        tmp.append(np.random.randint(0, self.Lattice))
                    x_site = tmp[1]
                    y_site = tmp[4]
                    top = (x_site - 1) % self.Lattice
                    bottom = (x_site + 1) % self.Lattice
                    left = (y_site - 1) % self.Lattice
                    right = (y_site + 1) % self.Lattice
                    delta_E = 2.0 * spin[x_site,y_site] * (self.J * (spin[top,y_site] + spin[bottom,y_site] + spin[x_site,left] + spin[x_site,right]) +  self.H )
                    if delta_E <= 0:
                        spin[x_site,y_site] = -spin[x_site,y_site]
                    else:
                        if np.exp((-delta_E)/t) > np.random.random():
                            spin[x_site,y_site] = -spin[x_site,y_site]
                                
            M = np.abs(np.mean(spin))
            M_2 = np.square(M)
            mag_2.append(M_2)
            mag.append(M)
            
        result = np.mean(mag[self.relax:])
        print("result is:", result)
        test = [] 

        #print("result is:", spin)
        result_2 = np.mean(mag_2[self.relax:]) 
        result_C = (result_2 - result ** 2) / T # 求磁比热
        test.append(t)
        test.append(result)
        test.append(result_C)
        
        np.save("M_T/{:d}.npy".format(T), test)
             
            #f_M.write(str(T)+'  '+str(result)+"\n") # 
            #f_C.write(str(T)+'  '+str(result_C)+"\n")
        
        #f_M.close()
        #f_C.close()
    def multi_run(self):
        
        start = time.time()
        
        ing_argv = [] 
        for T in range(8, 40, 1):
            ing_argv.append(T)
        
        with Pool(16) as p:
            p.map(self.main_loop, ing_argv)
        
        end = time.time() 
        
        print("run time is:", (end - start))
    
    def single_run(self):
        start = time.time()
        for T in range(3, 30, 1):
            self.main_loop(T)
        
        end = time.time() 
        print("run time is:", (end - start))

if __name__=="__main__":
    MC_model = MonteCarlo(100, 2000, 1, 1, 0, 10)
    MC_model.single_run()
