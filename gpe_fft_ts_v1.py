
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
# _______________________________________________________________________________

# Import libraries and general definitions
# -------------------------------------------------------------------------------

# In[1]:

import numpy as np
#import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
from scipy.integrate import odeint
#import numpy.linalg as lin
from gpe_fft_utilities import * # local folder utilities
from wave_functions import * # local folder initial states
# comment next line to export as a python shell 
#get_ipython().magic('matplotlib inline')
pi=np.pi


# Data block
# --------------------------------------------------------------------------------

# In[2]:

Zmax = 10      # Grid half length
Npoint = 128     # Number of grid points
Nparticle = 20   # Number of particles
a_s = 0.5        # scattering length 
whoz = 1.0       # harmonic oscilator angular frequency
Omega = 0* pi/(2*Zmax)    # reference frame velocity
Ntime_fin = 10000       # total number of time steps
Ntime_out = 100       # number of time steps for intermediate outputs
Dtr = 1.0e-3           # real time step
Dti = 1.0e-3           # imaginary time step
#
# print evolution data:
#
print("Initial data:")
print(" Number of particles = %g"%(Nparticle))
print(" Harmonic oscillator angular frequency = %g"%(whoz))
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
gint = 2*a_s*Nparticle*NormWF   # interaction (nonlinear-term) strength
Dz = 2*Zmax*1./Npoint              # length step size
Dk = pi/Zmax            # momentum step size
Kmax = Dk*(Npoint//2)   # maximum momentum
Dt = Dtr-1j*Dti         # complex time
Ninter = Ntime_fin//Ntime_out # Number of outputs with the intermediate states
print(" Characteristic interaction energy = %g"%(gint))


# Grid definitions: physical and momentum space
# ---------------------------------------------------------------------------------------

# In[4]:
print("dz = %g"%(Dz))

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
T_K = np.exp(-1j*0.5*Dt*Ekin_K)    # time evolution operator in K space (for second order accuracy)
# print("Ekin: "); print(Ekin_K)
#
# Potential energy in R space:
# Harmonic oscillator with angular frequency whoz:
Vpot_R = 0.5*whoz**2*zp**2;  
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


# Choose initial state (wave function)
# ---------------------------------------------------------------------------------------

# In[7]:

# initial wf: gaussian(x,n,x0,s0,w,v) 
# Gaussian centered at x0=0, width=w=1, velocity=v= 0, without nodes (s0=0)
c0=normaliza(gaussian(zp,Npoint,0,0,1,0)); # wf at t=0
#
# initial wf: thomas_fermi(x,n,x0,s0,gN,Ve) 
# Thomas Fermi for harmonic oscillator without nodes (s0=0)
#c0=normaliza(thomas_fermi(zp,Npoint,0,0,gint/NormWF,Vpot_R)); # wf at t=0


#  Evolve in time the initial state
#  ---------------------------------------------------------------------------------------

# In[8]:

# parameters for the time evolution
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
for i in range(1, Ntime_fin+1): # time evolution cicle  
    t += Dt.real
    psi=ifft(T_K*c)*Npoint
    c=T_K*fft(T_R_psi(t0,Dt,psi))/Npoint
    c = normaliza(c); # check norm in the wf
    if(not(i%Ntime_out)):
        j+=1
        tevol[j] = t
        energy_cicle[j,:] = Energy(c)
print("         final   = %g %g %g %g %g"%(Energy(c))) # check energies
print("Energy change at last step  = %g"%(energy_cicle[Ninter,0]-energy_cicle[Ninter-1,0]))


# Plot convergence during the evolution in the average energy per particle

# In[9]:

plot_convergence(tevol,energy_cicle[:,0],Ninter) # convergence in the average energy per particle


# Plot the final density (or wave function)

# In[10]:

cc = ifft(c)*Npoint*NormWF**0.5         # FFT from K3 to R3 and include the wf norm
psi = changeFFTposition(cc,Npoint,0)    # psi is the final wave function
plot_density(z,psi,Zmax,t)    
plot_phase(z,psi,Zmax,t)  
plot_real_imag(z,psi,Zmax,t)
cc0 = ifft(c0)*Npoint*NormWF**0.5         # FFT from K3 to R3 and include the wf norm
psi0 = changeFFTposition(cc0,Npoint,0)
plot_density(z,psi0,Zmax,0)
plot_real_imag(z,psi0,Zmax,0)

