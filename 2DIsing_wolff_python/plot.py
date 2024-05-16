import numpy as np 
import matplotlib.pyplot as plt
import math  
from scipy.optimize import curve_fit
import os, glob 

H = [0.0, 0.05, 0.1]

def load_npy(root):
    
    

    path = []
    
    for i in range(len(H)):
        path.append(os.path.join(root + "/M_T_H_{}".format(H[i]), "M_T.txt"))

    print("path is:", path)
    M = []
    T = []

    for i in range(len(path)):
        M_h, T_h = [], []
        data = np.loadtxt(path[i])
        for i in range(len(data)):
            T_h.append(data[i][0])
            M_h.append(data[i][1])
        M.append(M_h)
        T.append(T_h)
    
    print("M is:", M)
    
    return T, M

def load_dat(data_path): 
    
    T, M, C = [], [], [] 
    
    for line in open(data_path, "r"):
        tmp = line.split( )
        T.append(float(tmp[0]))
        M.append(float(tmp[1]))
    
    return T, M
    

def plot(T, M):
    
    for i in range(len(T)): 
        plt.plot(T[i],M[i],"o-", label="H="+str(H[i]))
    
    plt.legend(loc="upper right", prop = {'family': "Times New Roman", "weight":"normal", "size":16,}, frameon=False)
    
    plt.show()

if __name__=="__main__":
    path = "D:/Monte-Carlo-Algorithm-master/Monte-Carlo-Algorithm-master/2DIsing_wolff_python"
    T, M = load_npy(path)
    plot(T, M)