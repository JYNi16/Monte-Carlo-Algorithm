## Author jinyang ni 
## Monte-Carlo simulation for 2D-Ising model 

import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool 
from numba import jit 
import os

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
        t = T
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
            
        M_ave = np.mean(mag[self.relax:])
        print("M_ave is:", M_ave)

        #print("result is:", spin)
        M_2 = np.mean(mag_2[self.relax:]) 
        C_v = (M_2 - M_ave ** 2) / T # 求磁比热
        
        return M_ave
        #np.save("M_T/{:d}.npy".format(T), test)
             
            #f_M.write(str(T)+'  '+str(result)+"\n") # 
            #f_C.write(str(T)+'  '+str(result_C)+"\n")
        
        #f_M.close()
        #f_C.close()
    
    def save_data(self, data):
        save_path = "M_T"
        if os.path.exists(save_path):
            print("save path {} exist".format(save_path))
        else:
            print("save path {} not exist".format(save_path))
            os.makedirs(save_path)
            print("now makefir the save_path")

        np.savetxt(save_path + "/M_T.txt", data)
    
    def multi_run(self):
        
        start = time.time()
        
        M = []
        ing_argv = [] 
        T = [round(t,2) for t in np.linspace(0.8, 4, 33)]
        for i in range(len(T)):
            ing_argv.append(T[i])
        
        with Pool(8) as p:
            M.append(p.map(self.main_loop, ing_argv))
        
        print(M[0])
        
        data = [] 
        for i in range(len(T)):
            data.append([T[i], M[0][i]])
        
        self.save_data(data)
    
    def single_run(self):
        start = time.time()
        T = [round(t,2) for t in np.linspace(0.8, 4, 33)]
        M = []
        for i in range(len(T)):
            M.append(self.main_loop(T[i]))
        
        data = [] 
        print("M is:", M)
        for i in range(len(T)):
            data.append([T[i], M[i]])
        
        self.save_data(data)

if __name__=="__main__":
    MC_model = MonteCarlo(100, 2000, 1, 1, 0, 10)
    MC_model.single_run()
