import numpy as np 
import matplotlib.pyplot as plt
import math  
from scipy.optimize import curve_fit
import os, glob 

def load_MT(root):
    
    path_list = []
    path_list.extend(glob.glob(os.path.join(root, "*.npy")))
    
    M, T = [], []
    path_list.sort(key=len)
    
    for i in range(len(path_list)):
        data = np.load(path_list[i])
        print(data)
        M.append(data[1])
        T.append(data[0])
    
    return T, M 

def plot(T, M):
    plt.figure(figsize = (7,5))
    #p1 = plt.subplot(121)
    plt.plot(T,M,'o-', color = "#FF0000", markersize = 8, linewidth =1, label = r'$M$')
    #p1 = plt.subplot(122)
    #p1.plot(T,C,'o-', color = "#000000", markersize = 8, linewidth =1, label = r'$M$')
    plt.show() 


if __name__=="__main__":
    T, M = load_MT("E:/MonteCarlo/2D-Ising-Model/python/M_T")
    plot(T, M)
    
        
    
    


