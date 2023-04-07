from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

def q1(x, t):#defining the differential equation for the first question : 4x''+2500x=0
    return (x[1],-2500*x[0]/4-0*x[1]/4)

def q2(x, t):#defining the differential equation for the second question : 4x''+100x'+2500x=0
    return (x[1],-2500*x[0]/4-100*x[1]/4)

def q3(x, t):#defining the differential equation for the third question : 4x''+200x'+2500x=0
    return (x[1],-2500*x[0]/4-200*x[1]/4)

def q4(x, t):#defining the differential equation for the first question : 4x''+400x'+2500x=0
    return (x[1],-2500*x[0]/4-400*x[1]/4)


x0 = [0.15, 20]
ts=np.linspace(0, 5, 1000)
us = odeint( q1, x0, ts)
ys=us[:,0]
amplitude=ys.max()
plt.figure("1.1")
plt.plot(ts, ys,'-')
#plt.plot(ts,ys,'ro')
plt.xlabel('t')
plt.ylabel('x')
plt.show()
print("The amplitude of the vibration is",amplitude)

