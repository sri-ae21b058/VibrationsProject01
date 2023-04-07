from matplotlib import pyplot as plt
import numpy as np
import question1 as q1
def steady_state_prop(F0, k, m, c,force_freq):
    # F0 = initial force, k = spring constant, m = mass, c = damping coefficient, t = time
    # returns the steady state amplitude and phase angle
    X_steady=F0/np.sqrt((k-m*force_freq**2)**2+(c*force_freq)**2)#steady state amplitude
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
    roots = q1.charac_roots(m,k,c)
    #we use code from question1.py to find the roots of the characteristic equation of the homogeneous equation : mx''+cx'+kx=0
    [X_steady,phi_steady]=steady_state_prop(F0, k, m, c,force_freq)
    if roots[0].imag != 0 and roots[1].imag != 0:#complex roots underdamped case
        '''we use the solution of the homogeneous equation x(t)=e^(lambda1t)(X_0cos(omega_d*t+phi_0))'''
        X_0=np.sqrt((x0-X_steady*np.cos(phi_steady))**2+((v0-X_steady*np.sin(phi_steady)*force_freq+roots[0].real*(x0-X_steady*np.cos(phi_steady)))/roots[0].imag)**2)
        phi_0=np.arctan((v0-X_steady*np.sin(phi_steady)*force_freq+roots[0].real*(x0-X_steady*np.cos(phi_steady)))/(roots[0].imag*(x0-X_steady*np.cos(phi_steady))))
        return X_0*np.exp(roots[0].real*t)*np.cos(roots[0].imag*t-phi_0)+X_steady*np.cos(force_freq*t+phi_steady)
    else:
        return 0
    
def plot_soln(F0, k, m, c,force_freq,t,x0,v0,fignum):
    # F0 = initial force, k = spring constant, m = mass, c = damping coefficient, t = time
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
    return




