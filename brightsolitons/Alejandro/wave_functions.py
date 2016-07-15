
# coding: utf-8

# In[1]:

import numpy as np
from scipy.fftpack import fft  #, ifft
import scipy.optimize as opt
from sys import exit


# Wave functions definitions
# -------------------------------------------------------------------------------------

# In[2]:

def gaussian(x,n,x0,s0,w,v): # Gaussian wf in  K3
    # x  = grid points
    # n  = number of grid points
    # x0 = Gaussian position
    # s0 = switch for a node: =1 imprint node, =0 no node     
    # w =  Gaussian width
    # v =  wavepacket speed         
    #
    print('Initial state: Gaussian')
    fx = np.exp(1j*v*x)*np.pi**0.25*np.exp(-0.5*((x-x0)/w)**2); # define the Gaussian in R3
    if (s0):
        fx *= x  # Equivalent to a first excited state in a harmonic trap
    return fft(fx)/n   # FFT to K3


# In[ ]:

def thomas_fermi(x,n,x0,s0,gN,Ve): # Thomas-Fermi (with/out node) ansatz for harmonic oscillator in K3
    # It assumes the harmonic oscillator frequency as 1,
    # although this function should be updated for generic frequency and intended for generic potentials
    # x  = grid points
    # n  = number of grid points
    # Ve =  external-potential vector 
    # gN = interaction strength times number of particles 
    # x0 = node position    
    # s0 = switch for imprinting a (hyperbolic tangent) node: =1 imprint node, =0 no node    
    #
    print('Initial state: Thomas Fermi (TF) ansatz')
    R=(1.5*gN)**(1.0/3.0)
    muTF=0.5*R**2
    print('  TF chemical potential = ', muTF,', TF radius =',R)
    fx = (muTF-Ve)/gN # fx is now the TF density (including the negative values)
    for i in range(n):
        if ( fx[i] > 0.0 ):
            fx[i] =fx[i]**0.5
        else: 
            fx[i] = 0.0
    if ( s0 ):  # imprint a node
        fx *= np.tanh((x-x0)*(muTF-Ve)**0.5)  # the local healing length is \hbar/(m (mu-V(x))) 
    return fft(fx)/n   # FFT to K3


# In[1]:

def dark_soliton(x,n,x0,gn,v): # tanh(x/\xi) wf
    # x  = grid points
    # n  = number of grid points
    # x0 = soliton position
    # gn = (positive) interaction energy or chemical potential  
    # v =  soliton speed   
    #
    if (gn <= 0):
        exit(" error: interaction must be positive...")    
    print('Initial state: dark soliton: hyperbolic tangent')
    c = (gn)**0.5 # speed of sound
    xi = 1.0/c  
    print('  healing length = ',xi)    
    b = v/c
    print('  velocity / speed of sound = ',b )  
    d= (1-b**2)**0.5  # velocity dependent width
    fx =1j*b+d*np.tanh(d*(x-x0)/xi)  # (complex) wave function in R3
    return fft(fx)/n    # FFT to K3


# In[ ]:

def bright_soliton(x,n,x0,gn,v): # 1/cosh(x/\xi) wf
    # x  = grid points
    # n  = number of grid points
    # x0 = soliton position
    # gn = (negative) interaction energy or double of the chemical potential  
    # v =  soliton speed     
    #
    if (gn >= 0):
        exit(" error: interaction must be negative...")
    print('Initial state: bright soliton: hyperbolic secant')
    psi0 = 0.5*np.sqrt(abs(gn))
    xi = 1.0/(2*abs(gn)*psi0**2)**0.5
    print('healing length=',xi)
    print('  velocity = ',v )
    fx = np.exp(1j*v*x)*psi0/np.cosh((x-x0)/(np.sqrt(2)*xi))   # wave function in R3
    return fft(fx)/n     # FFT to K3

