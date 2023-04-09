from matplotlib import pyplot as plt
import numpy as np
def charac_roots(m,k,c):
    #m = mass, k = spring constant, c = damping coefficient
    #returns the roots of the characteristic equation: lambda^2 + 2*lambda*c/m + k/m = 0
    roots = np.zeros(2,dtype=complex)
    roots[1] = (-c + np.emath.sqrt(c**2 - 4*m*k))/(2*m) #np.emath.sqrt gives the complex square root
    roots[0] = (-c - np.emath.sqrt(c**2 - 4*m*k))/(2*m)
    return roots

def sol_eqn(x0,v0,t,m,k,c):
    #x0 = initial position, v0 = initial velocity, t = time, m = mass, k = spring constant, c = damping coefficient
    #returns the displacement and velocity of the mass at time t
    roots = charac_roots(m,k,c)
    if roots[0].imag != 0 and roots[1].imag != 0:#complex roots case underdamped
        A = np.array([[1,1],[roots[0],roots[1]]],dtype=complex)
        B = np.array([x0,v0])
        const=np.linalg.inv(A).dot(B) #finding the constants of the solution using inverse matrix of A and multiplying by B
        return np.exp(roots[0].real*t)*(const[0]*np.cos(roots[0].imag*t)+const[1]*np.sin(roots[0].imag*t))#x(t)=e^(lambda1t)(C1*cos(omega1t)+C2*sin(omega1t))
    elif roots[0].imag == 0 and roots[1].imag == 0: 
        if roots[0].real == roots[1].real:#single root case Critically damped
            A=np.array([[1,0],[roots[0],roots[0]+1]])
            B=np.array([x0,v0])
            const=np.linalg.inv(A).dot(B) #finding the constants of the solution
            return const[0]*np.exp(roots[0].real*t)+const[1]*t*np.exp(roots[0].real*t)#x(t)=e^(lambda1t)(C1+C2t)
        else:# overdamped case
            A=np.array([[1,1],[roots[0],roots[1]]])
            B=np.array([x0,v0])
            const=np.linalg.inv(A).dot(B) #finding the constants of the solution
            return const[0]*np.exp(roots[0]*t)+const[1]*np.exp(roots[1]*t)#x(t)=C1*e^(lambda1t)+C2*e^(lambda2t)


def plot_sol(x0,v0,t,m,k,c,fignum): #plotting the solution, we only need the real part of the roots
    plt.figure(fignum)
    plt.clf()
    plt.plot(t,sol_eqn(x0,v0,t,m,k,c))
    plt.title('Damping coefficient = '+str(c))
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.grid()
    plt.savefig('question '+fignum+'.png')
    plt.show()

x0=0.15#given initial conditions
v0=-20
t=np.linspace(0,2,1000)
m=4#given parameters of mass, spring constant and damping coefficients
k=2500
c=[0,100,200,400]
nfig=['1.1','1.2','1.3','1.4']
for i in range(len(c)):
    plot_sol(x0,v0,t,m,k,c[i],nfig[i])
