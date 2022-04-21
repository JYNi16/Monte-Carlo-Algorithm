import numpy as np 
import matplotlib.pyplot as plt
import math  
from scipy.optimize import curve_fit

p = math.pi
T = [] 
M = [] 
C = []
i = 0 

with open ("M.dat",'r') as f1:
    while( i < 50):
        i = i + 1
        value = f1.readline()
        v = [float(s) for s in value.split()]
        T.append(v[0])
        M.append(v[1])
        C.append(v[2])


plt.figure(figsize = (16,5))
p1 = plt.subplot(121)
p1.plot(T,M,'o-', color = "#FF0000", markersize = 8, linewidth =1, label = r'$M$')
p1 = plt.subplot(122)
p1.plot(T,C,'o-', color = "#000000", markersize = 8, linewidth =1, label = r'$M$')
plt.show() 