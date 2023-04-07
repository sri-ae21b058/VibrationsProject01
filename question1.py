from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

def diff_eqn(x, t,k,c,m):#defining the differential equation for the first question : 4x''+2500x=0
    return (x[1],-k*x[0]/m-c*x[1]/m)


x0 = [0.15,-20]
ts=np.linspace(0, 2, 1000)
us = odeint( diff_eqn, x0, ts,args=(2500,0,4))
ys=us[:,0]
amplitude=ys.max()
plt.figure("1.1")
plt.plot(ts, ys,'-')
#plt.plot(ts,ys,'ro')
plt.xlabel('t')
plt.ylabel('x')
plt.show()
print("The amplitude of the vibration is",amplitude)
#for question 1.2
us = odeint( diff_eqn, x0, ts,args=(2500,100,4))
ys=us[:,0]
amplitude=ys.max()
plt.figure("1.2")
plt.plot(ts, ys,'-')
#plt.plot(ts,ys,'ro')
plt.xlabel('t')
plt.ylabel('x')
plt.show()
print("The amplitude of the vibration is",amplitude)
#for question 1.3
us = odeint( diff_eqn, x0, ts,args=(2500,200,4))
ys=us[:,0]
amplitude=ys.max()
plt.figure("1.3")
plt.plot(ts, ys,'-')
#plt.plot(ts,ys,'ro')
plt.xlabel('t')
plt.ylabel('x')
plt.show()
print("The amplitude of the vibration is",amplitude)
#for question 1.4
us = odeint( diff_eqn, x0, ts,args=(2500,400,4))
ys=us[:,0]
amplitude=ys.max()
plt.figure("1.4")
plt.plot(ts, ys,'-')
#plt.plot(ts,ys,'ro')
plt.xlabel('t')
plt.ylabel('x')
plt.show()
print("The amplitude of the vibration is",amplitude)


