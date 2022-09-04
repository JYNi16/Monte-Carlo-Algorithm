from xymodel_mc import *
import numpy as np

L = 50
beta = 1
J = 5
steps = 1000000
V = []
M2 = []
C = []

def single_run():
    for T in range(4, 25):
        t = T/10
        sim = XYMetropolis((L, L),beta=1/t,J=J,random_state=5,initial_state='hot')
        sim.simulate(steps)
        V.append(sim.Vdensity)
        M2.append(sim.M2)
        C.append(sim.C)

if __name__=="__mian__":
    single_run()