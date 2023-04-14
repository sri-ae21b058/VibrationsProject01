from matplotlib import pyplot as plt
import numpy as np
x0=[[1,2],[1,0]]#initial displacements
v0=[[1,-2],[-1,0]]#initial velocities
m=[2,1]#masses
k=[1000,500]#spring constants

def get_modal_freq(m,k):
    #returns the modal frequencies for the system of equations : m[0]*x''[0]+(k[0]+k[1])*x[0]-k[1]*x[1]=0, m[1]*x''[1]+(k[1])*x[1]-k[1]*x[0]=0
    f_modal=np.zeros(2)
    M=np.array([[m[0],0],[0,m[1]]])
    K=np.array([[k[0]+k[1],-k[1]],[-k[1],k[1]]])
    f_modal=np.sqrt(np.linalg.eigvals(np.linalg.inv(M).dot(K)))
    return f_modal

def get_modal_shapes(x0,v0,m,k,t,f_modal):
    #x0 = initial position, v0 = initial velocity, t = time, m = mass, k = spring constant, f_modal = modal frequencies
    #returns the modal shapes
    r=np.zeros(2)
    for i in range(2):
        r[i]=((k[1]+k[0])-m[0]*f_modal[i]**2)/k[1]
    A=np.array([[1,1],[r[0],r[1]]])
    B=np.array([x0[0],x0[1]])
    Xcos=np.linalg.inv(A).dot(B)
    A=np.array([[-f_modal[0],-f_modal[1]],[-r[0]*f_modal[0],-r[1]*f_modal[1]]])
    B=np.array([v0[0],v0[1]])
    Xsin=np.linalg.inv(A).dot(B)
    X=np.array([np.sqrt(Xcos[0]**2+Xsin[0]**2), np.sqrt(Xcos[1]**2+Xsin[1]**2)])
    if Xcos[0]==0: #in case of division by zero
        phi=np.array([np.pi/2, np.arctan(Xsin[1]/Xcos[1])])
    elif Xcos[1]==0:
        phi=np.array([np.arctan(Xsin[0]/Xcos[0]), np.pi/2])
    else:
        phi=np.array([np.arctan(Xsin[0]/Xcos[0]), np.arctan(Xsin[1]/Xcos[1])])
    x=np.zeros(2,dtype=object)
    x[0]=X[0]*np.cos(f_modal[0]*t+phi[0])+X[1]*np.cos(f_modal[1]*t+phi[1])
    x[1]=r[0]*X[0]*np.cos(f_modal[0]*t+phi[0])+r[1]*X[1]*np.cos(f_modal[1]*t+phi[1])
    return x
def plot_modal_shapes(x0,v0,m,k,t,f_modal,fignum):
    #plots the modal shapes
    x=get_modal_shapes(x0,v0,m,k,t,f_modal)
    plt.figure(fignum+' Displacement vs Time')
    plt.plot(t,x[0],label='mass 1')
    plt.plot(t,x[1],label='mass 2')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title(fignum+': Displacement vs Time')
    plt.grid()
    plt.savefig(fignum+' Displacement vs Time.png')
    plt.show()
    
f_modal=get_modal_freq(m,k)#the modal frequencies
t_modalmax=(2*np.pi/f_modal).max()#the modal period
t=np.linspace(0,t_modalmax,1000)#time domain
plot_modal_shapes(x0[0],v0[0],m,k,t,f_modal,'3.1')#plotting the modal shapes for x0=[1,2] and v0=[1,-2]
plot_modal_shapes(x0[1],v0[1],m,k,t,f_modal,'3.2')#plotting the modal shapes for x0=[1,0] and v0=[-1,0]
