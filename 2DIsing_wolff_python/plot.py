import numpy as np 
import matplotlib.pyplot as plt
import math  
from scipy.optimize import curve_fit
import os, glob 

def load_npy(root):
    
    path_list_0 = []
    path_list_1 = []
    path_list_2 = []
    path_list_3 = []
    path_list_4 = []
    path_list_5 = []

    path = []

    path_list_0.extend(glob.glob(os.path.join(root + "/M_T_0.0", "*.npy")))
    path_list_1.extend(glob.glob(os.path.join(root + "/M_T_0.01", "*.npy")))
    path_list_2.extend(glob.glob(os.path.join(root + "/M_T_0.02", "*.npy")))
    path_list_3.extend(glob.glob(os.path.join(root + "/M_T_0.03", "*.npy")))
    path_list_4.extend(glob.glob(os.path.join(root + "/M_T_0.04", "*.npy")))
    path_list_5.extend(glob.glob(os.path.join(root + "/M_T_0.05", "*.npy")))

    path = [path_list_0, path_list_1, path_list_2, path_list_3, path_list_4, path_list_5]

    M = [[], [], [], [], [], []]
    T = [[], [], [], [], [], []]

    for i in range(6):
        for j in range(len(path_list_0)):
            data = np.load(path[i][j])
            T[i].append(data[0])
            M[i].append(data[1])

    return T, M

def load_dat(data_path): 
    
    T, M, C = [], [], [] 
    
    for line in open(data_path, "r"):
        tmp = line.split( )
        T.append(float(tmp[0]))
        M.append(float(tmp[1]))
    
    return T, M
    

def plot(T, M):

    plt.plot(T[0],M[0],'o-', color = "b", markersize=6, linewidth =1, label=0.0)
    plt.plot(T[1], M[1], 'o-', color="g", markersize=6, linewidth=1, label=0.1)
    plt.plot(T[2], M[2], 'o-', color="r", markersize=6, linewidth=1, label=0.2)
    plt.plot(T[3], M[3], 'o-', color="c", markersize=6, linewidth=1, label=0.3)
    plt.plot(T[4], M[4], 'o-', color="m", markersize=6, linewidth=1, label=0.4)
    plt.plot(T[5], M[5], 'o-', color="y", markersize=6, linewidth=1, label=0.5)

    plt.show()

if __name__=="__main__":
    path = "E:/MonteCarlo/2D-Ising-Model/Wolff_cluster/2DIsing_wolff_python"
    T, M = load_npy(path)
    plot(T, M)