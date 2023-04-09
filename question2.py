from matplotlib import pyplot as plt
import numpy as np
def charac_roots(m,k,c):
    #m = mass, k = spring constant, c = damping coefficient
    #returns the roots of the characteristic equation: lambda^2 + 2*lambda*c/m + k/m = 0
    roots = np.zeros(2,dtype=complex)
    roots[1] = (-c + np.emath.sqrt(c**2 - 4*m*k))/(2*m)
    roots[0] = (-c - np.emath.sqrt(c**2 - 4*m*k))/(2*m)
    return roots
def steady_state_prop(F0, k, m, c,force_freq):
    # F0 = initial force, k = spring constant, m = mass, c = damping coefficient, t = time
    # returns the steady state amplitude and phase angle
    X_steady=F0/np.sqrt((k-m*force_freq**2)**2+(c*force_freq)**2)#steady state amplitude
    if k-m*force_freq**2==0:#if the denominator is zero, we use the limit of the ratio as the phase angle
        phi_steady=np.pi/2
    else:
        phi_steady=np.arctan((c*force_freq)/(k-m*force_freq**2))#steady state phase angle
    return X_steady,phi_steady

def steady_state_soln(F0, k, m, c,force_freq,t):
    # F0 = initial force, k = spring constant, m = mass, c = damping coefficient, t = time
    # returns the steady state solution at time t
    [X_steady,phi_steady]=steady_state_prop(F0, k, m, c,force_freq)
    return X_steady*np.cos(force_freq*t+phi_steady)

def total_soln(F0, k, m, c,force_freq,t,x0,v0):
    # F0 = initial force, k = spring constant, m = mass, c = damping coefficient, t = time
    # returns the total force at time t
    roots = charac_roots(m,k,c)
    #we use code from question1.py to find the roots of the characteristic equation of the homogeneous equation : mx''+cx'+kx=0
    [X_steady,phi_steady]=steady_state_prop(F0, k, m, c,force_freq)
    if roots[0].imag != 0 and roots[1].imag != 0:#complex roots underdamped case
        '''we use the solution of the homogeneous equation x(t)=e^(lambda1t)(X_0cos(omega_d*t+phi_0))'''
        X_0=np.sqrt((x0-X_steady*np.cos(phi_steady))**2+((v0-X_steady*np.sin(phi_steady)*force_freq+roots[0].real*(x0-X_steady*np.cos(phi_steady)))/roots[0].imag)**2)
        if (roots[0].imag*(x0-X_steady*np.cos(phi_steady)))==0:#if the denominator is zero, we use the limit of the ratio as the phase angle
            phi_0=np.pi/2
            return X_0*np.exp(roots[0].real*t)*np.sin(roots[0].imag*t)+X_steady*np.cos(force_freq*t+phi_steady)
        else:
            phi_0=np.arctan((v0-X_steady*np.sin(phi_steady)*force_freq+roots[0].real*(x0-X_steady*np.cos(phi_steady)))/(roots[0].imag*(x0-X_steady*np.cos(phi_steady))))
            return X_0*np.exp(roots[0].real*t)*np.cos(roots[0].imag*t-phi_0)+X_steady*np.cos(force_freq*t+phi_steady)
    else:
        return 0
    
def plot_soln(F0, k, m, c,force_freq,t,x0,v0,fignum):
    # F0 = initial force, k = spring constant, m = mass, c = damping coefficient, t = time, fignum = figure number, x0 = initial displacement, v0 = initial velocity
    # plots the total force at time t
    plt.figure(fignum+'total solution')
    plt.plot(t,total_soln(F0, k, m, c,force_freq,t,x0,v0))
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title('Total solution of the forced damped harmonic vibration F(t)='+str(F0)+'cos('+str(force_freq)+'t)')
    plt.savefig(fignum+'total solution.png')
    plt.show()
    plt.figure(fignum+'steady state solution')
    plt.plot(t,steady_state_soln(F0, k, m, c,force_freq,t))
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.title('Steady state solution of the forced damped harmonic vibration F(t)='+str(F0)+'cos('+str(force_freq)+'t)')
    plt.savefig(fignum+'steady state solution.png')
    plt.show()

x0=0.1#initial displacement
v0=10#initial velocity
m=10#mass
k=4000#spring constant
c=40#damping coefficient
F0=100#Force amplitude
force_freq=[10,20]#frequency of the force
t=np.linspace(0,2,1000)#time
nfig=['2.1','2.2']
for i in range(2):
    plot_soln(F0, k, m, c,force_freq[i],t,x0,v0,nfig[i])





