from matplotlib import pyplot as plt
import numpy as np
from sympy import *
x0=[[1,2],[1,0]]#initial displacements
v0=[[1,-2],[-1,0]]#initial velocities
m=[1,2]#masses
k=[1000,500]#spring constants

def onlypositive(l):
    #returns only the positive roots
    l1=[]
    for i in l:
        if i>0:
            l1.append(i)
    return l1
def get_modal_freq(m,k):
    #returns the modal frequencies
    f=Symbol('f')
    A=Matrix([[(k[1]+k[0])-m[0]*f**2,-k[1]],[-k[1],k[1]-m[1]*f**2]])
    f_modal=onlypositive(solve(det(A),f))#solving the fourth order equation from det(A)=0 and taking only the positive roots gives us the two eigen values we need
    return f_modal

def get_modal_shapes(x0,v0,m,k,t,f_modal):
    #x0 = initial position, v0 = initial velocity, t = time, m = mass, k = spring constant, f_modal = modal frequencies
    #returns the modal shapes
    
    r=zeros(2)
    for i in range(2):
        r[i]=(k[1]+k[0])-m[0]*f_modal[i]**2/-k[1]
    A=Matrix([[1,1],[r[0],r[1]]])
    B=Matrix([x0[0],x0[1]])
    Xcos=A.inv()*B
    A=Matrix([[-f_modal[0],-f_modal[1]],[-r[0]*f_modal[0],-r[1]*f_modal[1]]])
    B=Matrix([v0[0],v0[1]])
    Xsin=A.inv()*B
    X=np.array([np.sqrt(Xcos[0]**2+Xsin[0]**2), np.sqrt(Xcos[1]**2+Xsin[1]**2)])
    phi=np.array([np.arctan(Xsin[0]/Xcos[0]), np.arctan(Xsin[1]/Xcos[1])])
    x=zeros(2)
    x[0]=X[0]*np.cos(f_modal[0]*t+phi[0])+X[1]*np.cos(f_modal[1]*t+phi[1])
    x[1]=r[0]*X[0]*np.cos(f_modal[0]*t+phi[0])+r[1]*X[1]*np.cos(f_modal[1]*t+phi[1])
    v=zeros(2)
    v[0]=-f_modal[0]*X[0]*np.sin(f_modal[0]*t+phi[0])-f_modal[1]*X[1]*np.sin(f_modal[1]*t+phi[1])
    v[1]=-f_modal[0]*r[0]*X[0]*np.sin(f_modal[0]*t+phi[0])-f_modal[1]*r[1]*X[1]*np.sin(f_modal[1]*t+phi[1])
    return x,v

def plot_modal_shapes(x0,v0,m,k,t,f_modal,fignum):
    #plots the modal shapes
    [x,v]=get_modal_shapes(x0,v0,m,k,t,f_modal)
    plt.figure(fignum)
    plt.subplot(211)
    plt.plot(t,x[0],label='mass 1')
    plt.plot(t,x[1],label='mass 2')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title(fignum+': Displacement vs Time')
    plt.subplot(212)
    plt.plot(t,v[0],label='mass 1')
    plt.plot(t,v[1],label='mass 2')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.title(fignum+'Velocity vs Time')
    plt.show()

f_modal=get_modal_freq(m,k)#the modal frequencies
t=np.linspace(0,0.1,1000)#time domain
plot_modal_shapes(x0[0],v0[0],m,k,t,f_modal,'3.1')#plotting the modal shapes for x0=[1,2] and v0=[1,-2]
plot_modal_shapes(x0[1],v0[1],m,k,t,f_modal,'3.2')#plotting the modal shapes for x0=[1,0] and v0=[-1,0]
