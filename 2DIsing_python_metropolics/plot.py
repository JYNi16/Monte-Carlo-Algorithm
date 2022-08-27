import numpy as np 
import matplotlib.pyplot as plt
import math  
from scipy.optimize import curve_fit
import os, glob 

def load_mpy(root):
    
    path_list = []
    path_list.extend(glob.glob(os.path.join(root, "*.npy")))
    
    M, T, C = [], [], []
    path_list.sort(key=len)
    
    for i in range(len(path_list)):
        data = np.load(path_list[i])
        print(data)
        M.append(data[1])
        T.append(data[0])
        C.append(data[2])
    
    return T, M, C

def load_dat(data_path): 
    
    T, M, C = [], [], [] 
    
    for line in open(data_path, "r"):
        tmp = line.split( )
        T.append(float(tmp[0]))
        M.append(float(tmp[1]))
        C.append(float(tmp[2]))
    
    return T, M, C
    

def plot(T, M, C):
    plt.figure(figsize = (16,5))
    p1 = plt.subplot(121)
    plt.plot(T,M,'o-', color = "#FF0000", markersize = 8, linewidth =1, label = r'$M$')
    p1 = plt.subplot(122)
    p1.plot(T,C,'o-', color = "#000000", markersize = 8, linewidth =1, label = r'$M$')
    #plt.show()
    plt.savefig("test.jpg", dpi=500)


if __name__=="__main__":
    T,M,C = load_dat("E:/MonteCarlo/2D-Ising-Model/C++/m1.dat")
    plot(T, M, C)
    
        
    
    


