import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from numba import jit
import time
import os

sweeps = 600
K = 1
J = 1
H = 0
Lattice = 30
relax = 200

def down(x, y):
    if y > Lattice - 1.5: return [x, 0]
    else: return [x, y + 1]

@jit
def wolff_cluster(spin, P_add):
    mag = []
    mag_2 = []

    for k in range(sweeps + relax):
        x = np.random.randint(0, Lattice)
        y = np.random.randint(0, Lattice)
        sign = spin[x, y]
        stack = [[x, y]]
        lable = [[1 for i in range(Lattice)] for j in range(Lattice)]
        lable[x][y] = 0

        while len(stack) > 0.5:

            # While stack is not empty, pop and flip a spin
            [x_site, y_site] = stack.pop()
            spin[x_site, y_site] = -sign

            # Append neighbor spins

            # Left neighbor
            if x_site < 0.5:
                [leftx, lefty] = [Lattice - 1, y_site]
            else:
                [leftx, lefty] = [x_site - 1, y_site]

            if spin[leftx, lefty] * sign > 0.5 and np.random.rand() < P_add:
                stack.append([leftx, lefty])
                lable[leftx][lefty] = 0

            # Right neighbor
            if x_site > Lattice - 1.5:
                [rightx, righty] = [0, y_site]
            else:
                [rightx, righty] = [x_site + 1, y_site]

            if spin[rightx, righty] * sign > 0.5 and np.random.rand() < P_add:
                stack.append([rightx, righty])
                lable[rightx][righty] = 0

            # Up neighbor
            if y_site < 0.5:
                [upx, upy] = [x_site, Lattice - 1]
            else:
                [upx, upy] = [x_site, y_site - 1]

            if spin[upx, upy] * sign > 0.5 and np.random.rand() < P_add:
                stack.append([upx, upy])
                lable[upx][upy] = 0

            # Down neighbor
            if y_site > Lattice - 1.5:
                [downx, downy] = [x_site, 0]
            else:
                [downx, downy] = [x_site, y_site + 1]

            if spin[downx, downy] * sign > 0.5 and np.random.rand() < P_add:
                stack.append([downx, downy])
                lable[downx][downy] = 0

        if k > relax:
            M = np.abs(np.mean(spin))
            M_2 = np.square(M)
            mag_2.append(M_2)
            mag.append(M)

    return mag, mag_2


def main_loop(T):
    t = T
    P_add = 1 - np.exp(-2 * J / t)
    spin = np.random.randint(-1, 1, (Lattice, Lattice))
    spin[spin == 0] = 1
    mag, mag_2 = wolff_cluster(spin, P_add)
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


if __name__ == "__main__":
    multi_run()