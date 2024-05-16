import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from numba import jit
import time
import os

sweeps = 500
K = 1
J = 1
H = 0.1
Lattice = 40
relax = 100

font = {'family': "Times New Roman", "weight":"normal", "size":16,}

def left(x, y):
    if x < 0.5: return [Lattice - 1, y]
    else: return [x - 1, y]

def right(x, y):
    if x > Lattice - 1.5: return [0, y]
    else: return [x + 1, y]

def up(x, y):
    if y < 0.5: return [x, Lattice - 1]
    else: return [x, y - 1]

def down(x, y):
    if y > Lattice - 1.5: return [x, 0]
    else: return [x, y + 1]

#@jit
def wolff_cluster(spin, t):
    mag = []
    mag_2 = []
    for k in range(sweeps + relax):
        x = np.random.randint(0, Lattice)
        y = np.random.randint(0, Lattice)
        sign = spin[x, y]
        P_add = 1 - np.exp(-(2*J + H) / t)
        stack = [[x, y]]
        lable = np.ones([Lattice, Lattice], int)
        lable[x, y] = 0

        while len(stack) > 0.5:

            # While stack is not empty, pop and flip a spin
            [x_site, y_site] = stack.pop()
            spin[x_site, y_site] = -sign

            # Append neighbor spins

            # Left neighbor

            [leftx, lefty] = left(x_site, y_site)

            if spin[leftx, lefty] * sign > 0.5  and np.random.rand() < P_add:
                stack.append([leftx, lefty])
                lable[leftx, lefty] = 0

            # Right neighbor

            [rightx, righty] = right(x_site, y_site)

            if spin[rightx, righty] * sign > 0.5 and np.random.rand() < P_add:
                stack.append([rightx, righty])
                lable[rightx, righty] = 0

            # Up neighbor

            [upx, upy] = up(x_site, y_site)

            if spin[upx, upy] * sign > 0.5 and np.random.rand() < P_add:
                stack.append([upx, upy])
                lable[upx, upy] = 0

            # Down neighbor

            [downx, downy] = down(x_site, y_site)

            if spin[downx, downy] * sign > 0.5 and  np.random.rand() < P_add:
                stack.append([downx, downy])
                lable[downx, downy] = 0

        if k > relax:
            M = np.abs(np.mean(spin))
            M_2 = np.square(M)
            mag_2.append(M_2)
            mag.append(M)

    return mag, mag_2


def main_loop(t):
    spin = np.ones([Lattice, Lattice], dtype=int)

    mag, mag_2 = wolff_cluster(spin, t)

    result = np.mean(mag)
    print("t is", t, "result is:", result)
    test = []

    # print("result is:", spin)
    result_2 = np.mean(mag_2)
    result_C = (result_2 - result ** 2) / t  # 求磁比热
    
    return result_2

def save_data(data):
    save_path = "M_T_H_{}".format(H)
    if os.path.exists(save_path):
        print("save path {} exist".format(save_path))
    else:
        print("save path {} not exist".format(save_path))
        os.makedirs(save_path)
        print("now makefir the save_path")

    np.savetxt(save_path + "/M_T.txt", data)


def single_run():
    start = time.time()
    for T in range(3, 30, 1):
        main_loop(T)

    end = time.time()
    print("run time is:", (end - start))

def plot(data):
    t = []
    kxy = []
    print(data)
    for i in range(len(data)):
        t.append(data[i][0])
        kxy.append(data[i][1])
    
    #print("t is:", t, "kxy is:", kxy)
    plt.plot(t, kxy, "o-")
    plt.show()

def multi_run():
    start = time.time()
    T = [round(t,2) for t in np.linspace(0.8, 4, 33)]
    M = []
    ing_argv = []
    for i in range(len(T)) :
        ing_argv.append(T[i])
    
    with Pool(8) as p:
        M.append(p.map(main_loop, ing_argv)) 

    print(M[0])
    
    data = [] 
    for i in range(len(T)):
        data.append([T[i], M[0][i]])
    
    save_data(data)
    plot(data)


if __name__ == "__main__":
    multi_run()