
# coding: utf-8

# ## FFT solver for 1D Gross-Pitaevski equation

# We look for the complex function $\psi(x)$ satisfying the GP equation
# 
# $ i\partial_t \psi = \frac{1}{2}(-i\partial_x - \Omega)^2\psi+ V(x)\psi + g|\psi|^2\psi $,
# 
# with periodic boundary conditions.
# 
# Integration: pseudospectral method with split (time evolution) operator; 
# that is evolving in real (R) or momentum (K) space according to the operators
# in the Hamiltonian, i.e.
# first we evaluate
# 
# $\hat{\psi}(x,\frac{t}{2})=\mathcal{F}^{-1}\left[\exp\left(-i \frac{\hbar^2 k^2}{2} \frac{t}{2}\right)\,\psi(k,0)\right] $
# 
# and later
# 
# $\psi(k,t) = \exp(-i \frac{\hbar^2 k^2}{2} \frac{t}{2})\, 
# \mathcal{F}\left[\exp\left(-i (V(x)+\lvert \hat{\psi}(x,\frac{t}{2})\rvert ^2)\, t \right)\,\hat{\psi}(x,\frac{t}{2}) \,
# \right]$
# 
# where $\cal{F}$ is the Fourier transform.
# _______________________________________________________________________________
# 
# 
#
# _______________________________________________________________________________

# Import libraries and general definitions
# -------------------------------------------------------------------------------

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
from scipy.integrate import odeint
#import numpy.linalg as lin
from gpe_fft_utilities import * # local folder utilities
from wave_functions import * # local folder initial states
# comment next line to export as a python shell 
#get_ipython().magic('matplotlib inline')
pi=np.pi

import matplotlib.animation as animation

# Data block
# --------------------------------------------------------------------------------

# In[2]:

Zmax = 2.0**7      # Grid half length
Npoint = 2**10     # Number of grid points - better if power of 2 
Nparticle = 20   # Number of particles
a_s = -0.01        # scattering length
x0 = -30.0            #soliton initial position
v0 = 0.0            #initial velocity
v = 10.0           #inital push
Omega = 0* pi/(2*Zmax)    # reference frame velocity
Ntime_fin = 20000       # total number of time steps
Ntime_out = 100       # number of time steps for intermediate outputs
Dtr = 1.0e-3*1.0          # real time step
Dti = 1.0e-3*100           # imaginary time step
#
# print evolution data:
#
print("Initial data:")
print(" Number of particles = %g"%(Nparticle))
print(" Domain half length = %g"%(Zmax))
print(" Number of grid points = %g"%(Npoint))
print(" Scattering length = %g"%(a_s))
print(" Total time of evolution = %g"%(Ntime_fin*Dtr))
print(" Real time step = %g"%(Dtr))
print(" Imaginary time = %g"%(Dti))
print(" Intermediate solutions = %g"%(Ntime_fin/Ntime_out-1))


# Derived quantities
# -------------------------------------------------------------------------------------

# In[3]:

NormWF = 1.0/(2*Zmax)           # Wave function (WF) norm
gn = 2*a_s*Nparticle
gint = gn*NormWF   # interaction (nonlinear-term) strength
xi = 1.0/(0.5*abs(gn)**2)**0.5         # healing length
Dz = 2*Zmax/Npoint              # length step size
Dk = pi/Zmax            # momentum step size
Kmax = Dk*(Npoint//2)   # maximum momentum
Dt = Dtr-1j*Dti         # complex time
Ninter = Ntime_fin//Ntime_out # Number of outputs with the intermediate states

print(" Characteristic interaction energy = %g"%(gn))

# Potential parameters
ww = 0.75*Zmax        # walls width at a %1 of the grid
wh = 0.0               #walls height
bx = 0.0*Zmax       #barrier position a %1 of the grid
bw = 1*xi          #barrier width a %1 of the healing leangth
bl = bx-bw*0.5      #barrier left wall
br = bx+bw*0.5      #barrier right wall
bh = 0.0            #initial barrier heigh  
a = 0.5             #final barrier factor a*Energy 
who = 0.0       # harmonic oscilator angular frequency          

# Grid definitions: physical and momentum space
# ---------------------------------------------------------------------------------------

# In[4]:

z = np.arange(-Zmax+Dz,Zmax+Dz,Dz) # physical (R-space) grid points in ascending order
# zp=[(0:Dz:Zmax) (-(Zmax-Dz):Dz:-Dz)]; 
zp = changeFFTposition(z,Npoint,1) # (R-space) grid points with FFT order
#print("grid points (K-order): "); print(zp)
#print(" R-order: "); print(z)
#
# kp=[(0:Dk:Kmax) (-(Kmax-Dk):Dk:-Dk)]; # grid points (K-space), FFT order
kp = np.arange(-Kmax+Dk,Kmax+Dk,Dk)
kp = changeFFTposition(kp,Npoint,1)
#print("momentum values: "); print(kp)


# Define operators
# ---------------------------------------------------------------------------------------

# In[5]:

Ekin_K = 0.5*(kp-Omega)**2            # Kinetic energy in K space
def T_K(t):                           # time evolution operator in K space (for second order accuracy)
    return np.exp(-1j*0.5*t*Ekin_K)    
# print("Ekin: "); print(Ekin_K)
#
# Potential energy in R space:
# Harmonic oscillator with angular frequency whoz:
# Vpot_R = 0.5*whoz**2*zp**2;  
# print("Vpot: "); print(Vpot_R)    
#

def potencial(bh):
    # using Heaviside step function:
    # f(x)=h/(1+exp(m(x-x0)))
    # Squared box with the grid size:
    m = 1000
    Vpot = 0.0
    Vpot = wh-wh/(1+np.exp(m*(z-ww)))+wh/(1+np.exp(m*(z+ww)))    
    Vbar = bh/(1+np.exp(m*(z-br)))-bh/(1+np.exp(m*(z-bl)))
    Vho = 0.5*who**2*z**2
    Vpot=Vpot+Vbar+Vho
    return Vpot
        # print("Vpot: "); print(Vpot_R)

    
    
# Main functions
# ---------------------------------------------------------------------------------------

# In[6]:

def Energy(c): # Energy (per particle) calculation
    global gint, Vpot_R, Ekin_K, Npoint
    ek = sum(Ekin_K*abs(c)**2)              # Kinetic energy in K
    psi = ifft(c)*Npoint;                   # wf FFT to R
    ep = sum(Vpot_R*abs(psi)**2)/Npoint;    # Potential energy
    ei = 0.5*gint*sum(abs(psi)**4)/Npoint;  # Interaction energy
    em =  ek+ep+ei;                         # average energy
    chem_pot = em+ei;                       # chemical potential
    return em, chem_pot, ek, ep, ei
#
def T_R_psi(t,Dt,psi): # Action of the time evolution operator over state c in R space
    global gint, Vpot_R
    # Includes the external potential and the interaction operators:
    #       T_R_psi = exp(-i Dt (Vpot_R+ gint|psi|^2) ) c    
    # psi is the wave function in R space
    # t is the time (which is not used for time independant Hamiltonians)
    # Dt is the complex time step 
    #
    return np.exp( -1j*Dt*(Vpot_R + gint*(abs(psi)**2)) )*psi    # return action on psi

def Integral(x):
    ileft = integral(x,Dz,z,-Zmax,bl)
    iin = integral(x,Dz,z,bl,br)
    iright = integral(x,Dz,z,br,Zmax)
    return ileft, iin, iright

# Choose initial state (wave function)
# ---------------------------------------------------------------------------------------

# In[7]:

# initial wf: gaussian(x,n,x0,s0,w,v) 
# Gaussian centered at x0=0, width=w=1, velocity=v= 0, without nodes (s0=0)
#c0=normaliza(gaussian(zp,Npoint,3,0,1,0)); # wf at t=0
#
# initial wf: thomas_fermi(x,n,x0,s0,gN,Ve) 
# Thomas Fermi for harmonic oscillator without nodes (s0=0)
#c0=normaliza(thomas_fermi(zp,Npoint,0,0,gint/NormWF,Vpot_R)); # wf at t=0
#
# initial wf: bright_soliton(x,n,x0,gn,v0): # 1/cosh(x/\xi)
# Bright Soliton 
c0=normaliza(bright_soliton(zp,Npoint,x0,gn,v0)) # wf at t=0

#  Evolve in time the initial state
#  ---------------------------------------------------------------------------------------

# In[8]:

# parameters for the imaginary time evolution
Vpot_R = changeFFTposition(potencial(bh),Npoint,1)     # prepare for fftpotencial(0)
t0=0.0 
tevol=np.empty([Ninter+1]) # time vector
energy_cicle=np.empty([Ninter+1,5]) # put the energies in a matrix
energy_cicle[0,:] = Energy(c0) # Energies at t=0
print("Energies:          Emed    mu    Ekin    Epot    Eint")
print("         initial = %g %g %g %g %g"%(Energy(c0)))
# print("$\psi(t=0)$: "); print(ct)  # check
c=c0 # initialize
tevol[0]=t0
j=0
t=0

#f3 = plt.figure()
cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function

plt.plot(z, abs(psi)**2, 'r-',label='$|\psi|^2$') # plot initial density
plt.plot(z, potencial(bh), 'g-') # plot the box

for i in range(1, Ntime_fin+1): # time evolution cicle  
    t += Dt.real
    psi=ifft(T_K(-1j*Dti)*c)*Npoint
    c=T_K(-1j*Dti)*fft(T_R_psi(t0,-1j*Dti,psi))/Npoint
    c = normaliza(c); # check norm in the wf
    if(not(i%Ntime_out)):
        j+=1
        tevol[j] = t
        energy_cicle[j,:] = Energy(c)
        
        if(not(j%10)):
            #prepare to plot
            cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
            psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function
            # plot features
            
            plt.title('Imaginary Time evolution'%(tevol[Ninter]),fontsize=15)
            plt.xlabel('$x/a_{ho}$',fontsize=15)
            #plt.ylabel('$\\psi\,(x)$',fontsize=15)
            plt.axis([-Zmax,Zmax,0, 0.3])
            plt.xticks(np.arange(-Zmax, Zmax+1,Zmax/2))
            plt.locator_params('y',nbins=3)
            #plt.plot(z, psi.real, 'r.',label='real$(\psi)$')
            #plt.plot(z, psi.imag, 'b--',label='imag$(\psi)$')
            plt.plot(z, abs(psi)**2, 'b--',label='$|\psi|^2$') # plot density
            #plt.legend(fontsize=15)
           #f3.show()        
            
print("         final   = %g %g %g %g %g"%(Energy(c))) # check energies
print("Energy change at last step  = %g"%(energy_cicle[Ninter,0]-energy_cicle[Ninter-1,0]))


# Plot convergence during the evolution in the average energy per particle
#  ---------------------------------------------------------------------------------------

# In[9]:

plot_convergence(tevol,energy_cicle[:,0],Ninter) # convergence in the average energy per particle

# Plot the final density (or wave function)
#  ---------------------------------------------------------------------------------------

# In[10]:

cc = ifft(c)*Npoint*NormWF**0.5         # FFT from K3 to R3 and include the wf norm
psi = changeFFTposition(cc,Npoint,0)    # psi is the final wave function
#plot_density(z,psi,Zmax,t)    
#plot_phase(z,psi,Zmax,t)  
#plot_real_imag(z,psi,Zmax,t)
cc0 = ifft(c0)*Npoint*NormWF**0.5         # FFT from K3 to R3 and include the wf norm
psi0 = changeFFTposition(cc0,Npoint,0)
#plot_density(z,psi0,Zmax,0)
#plot_real_imag(z,psi0,Zmax,0)

##write psi
#fn = open("densidad" , "w")
#format = "%g \t %g \n"
#for i in range(0, Npoint):
#    fn.write(format%(z[i], abs(psi[i])**2))
#fn.close()


#  Evolve in real time the eigenstate
#  ---------------------------------------------------------------------------------------

# In[11]:

# parameters for the real time evolution
psi_v = psi*np.exp(1j*v*z)        # giving a push to the eigenstate
c = fft(changeFFTposition(psi_v,Npoint,1))
c = normaliza(c); # check norm in the wf
c0=c
bh = Energy(c)[0]*a # barrier height
Vpot_R = changeFFTposition(potencial(bh),Npoint,1)     # prepare for fft     # Add the barrier. # Prepare for fft

t0=0.0
tevol=np.empty([Ninter+1]) # time vector
energy_cicle=np.empty([Ninter+1,5]) # put the energies in a matrix
energy_cicle[0,:] = Energy(c) # Energies at t=0
print("Energies:          Emed    mu    Ekin    Epot    Eint")
print("         initial = %g %g %g %g %g"%(Energy(c0)))
# print("$\psi(t=0)$: "); print(ct)  # check
tevol[0]=t0
j=0
t=0
f4 = plt.figure()
cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function

wave_function = np.empty([Ninter+1,3]) # put the % of the wf in a matrix
wave_function[0,:] = Integral(abs(psi)**2) # % at t=0

plt.plot(z, abs(psi)**2, 'r-',label='$|\psi|^2$') # plot initial density
plt.plot(z, potencial(bh), 'g-') # plot the box

for i in range(1, Ntime_fin+1): # time evolution cicle  
    t += Dt.real
    psi=ifft(T_K(Dtr)*c)*Npoint
    c=T_K(Dtr)*fft(T_R_psi(t0,Dtr,psi))/Npoint
    c = normaliza(c); # check norm in the wf
    if(not(i%Ntime_out)):
        j+=1
        tevol[j] = t
        energy_cicle[j,:] = Energy(c)
        cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
        psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function
        wave_function[j,:] = Integral(abs(psi)**2)  
        #print(wave_function[j,:])
        
        if(not(j%10)):
            #prepare to plot
        # comment next 2 lines if already done
#            cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
#            psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function
            # plot features
            
            plt.title('Real Time evolution'%(tevol[Ninter]),fontsize=15)
            plt.xlabel('$x/a_{ho}$',fontsize=15)
            #plt.ylabel('$\\psi\,(x)$',fontsize=15)
            plt.axis([-Zmax,Zmax,0, 0.3])
            plt.xticks(np.arange(-Zmax, Zmax+1,Zmax/2))
            plt.locator_params('y',nbins=3)
            #plt.plot(z, psi.real, 'r.',label='real$(\psi)$')
            #plt.plot(z, psi.imag, 'b--',label='imag$(\psi)$')
            plt.plot(z, abs(psi)**2, 'b--',label='$|\psi|^2$') # plot density
            #plt.legend(fontsize=15)
            f4.show()   
            
print("         final   = %g %g %g %g %g"%(Energy(c))) # check energies
print("Energy change at last step  = %g"%(energy_cicle[Ninter,0]-energy_cicle[Ninter-1,0]))

plot_convergence(tevol,energy_cicle[:,0],Ninter) # convergence in the average energy per particle
plot_wave_function(tevol, wave_function, Ninter) # % in time


#  Animation
#  ---------------------------------------------------------------------------------------

# In[12]:

#figure window
fig = plt.figure()
ax = plt.axes(xlim=(-Zmax, Zmax), ylim=(0, 0.3))
line, = ax.plot([], [], lw=2)
c=c0
cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function
#base frame
def init():
    line.set_data([], [])
    plt.plot(z, abs(psi)**2, 'r-',label='$|\psi_0|^2$') # plot initial density
    plt.plot(z, potencial(bh), 'g-') # plot the box
    return line,

# animation function.  This is called sequentially
def animate(i):
    global c
    psi=ifft(T_K(Dtr)*c)*Npoint
    c=T_K(Dtr)*fft(T_R_psi(t0,Dtr,psi))/Npoint
    c = normaliza(c); # check norm in the wf
    #prepare to plot
    cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
    psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function
    # plot features
    line.set_data(z, abs(psi)**2)
    return line,
    
#animation object
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=100, interval=20, blit=True)
fig.show()
plt.show()
