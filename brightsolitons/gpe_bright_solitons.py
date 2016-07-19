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
# ------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import os, glob
import matplotlib.animation as animation

from scipy.fftpack import fft, ifft
from scipy.integrate import odeint
from gpe_bs_utilities import *
from gpe_bs_parameters import *
pi=np.pi


# Print data block
# ------------------------------------------------------------------------------

print("\nInitial data:")
print(" Number of particles = %g"%(Nparticle))
print(" Harmonic oscillator angular frequency = %g"%(whoz))
print(" Domain half length = %g"%(Zmax))
print(" Number of grid points = %g"%(Npoint))
print(" Scattering length = %g"%(a_s))
print(" Total time of evolution = %g"%(Ntime_fin*Dtr))
print(" Real time step = %g"%(Dtr))
print(" Imaginary time = %g"%(Dti))
print(" Intermediate solutions = %g"%(Ntime_fin/Ntime_out-1))


# Print derived quantities and parameters of the initial wavefunction
# ------------------------------------------------------------------------------

print("\nInitial wavefunction parameters:")
print(" Characteristic interaction energy = %g"%(gint/NormWF))
print(" Healing length = %g \t %g" %(xi, 1.0 / ( np.abs(gn)**2 * 0.5   )**0.5))
print(" Position of the soliton = %g"%(x0))
if wall != 0.0:
    print(" Height of the walls = %g"%(wall))
if hb != 0.0:
    print(" Position of the potential barrier = %g" %(xb))
    print(" Height of the potential barrier = %g" %(hb))
    print(" Width of the potential barrier = %g" %(wb*2.0))
print(" Initial velocity of the soliton = %g \n" %(v/Zmax*0.5))


# Grid definitions: physical and momentum space; kinetic energy in K space
# ------------------------------------------------------------------------------

z = np.arange(-Zmax+Dz,Zmax+Dz,Dz)  # R-space grid points in ascending order
zp = changeFFTposition(z,Npoint,1)  # R-space grid points with FFT order
kp = np.arange(-Kmax+Dk,Kmax+Dk,Dk) # K-space grid points in ascending order
kp = changeFFTposition(kp,Npoint,1) # K-space grid points with FFT order
Ekin_K = 0.5*(kp-Omega)**2          # Kinetic energy in K space


# Initial state (bright soliton) and potential
# ------------------------------------------------------------------------------

t0=0.0

c0 = normaliza(bright_soliton(zp,Npoint,x0,gn),0);
c = c0 # initialize
psi = changeFFTposition(ifft(c)*Npoint*NormWF**0.5,Npoint,1)
psi0 = psi

# initial potential
Vpot_R = Vpot(0, z)

# saves the initial wavefunction (psi**2) and potential on a file
fv = open('initial1.dat', 'w')
format_v = "%.2f \t %.2f \t %.12g \n"
for i in range(0,Npoint-1):
    fv.write(format_v %(z[i],Vpot_R[i],np.abs(psi[i])**2))
fv.close()

# rearranges the potential with FFT order
Vpot_R = changeFFTposition(Vpot_R,Npoint,1) # fft order
V = changeFFTposition(Vpot_R,Npoint,0) # physical order



#  Evolve in time the initial state
# ------------------------------------------------------------------------------

# plots
plots = 0 # plots wavefunctions and intermediate states if 0

# checks evolution with imaginary time (comment next line to ignore it)
c = evolution(t0, -1j*Dti, z, c0, Vpot_R, V, Ekin_K, 1, plots)

# file writing (for real time evolution)
write_evolution = 1 # writes the wavefuntion at certain time steps if 0

# initial kick to the soliton (velocity v)
psi = psi *np.exp(+1j*v*(z-x0))
cc = changeFFTposition(psi,Npoint,1)
c = fft( cc / (Npoint*NormWF**0.5)); c = normaliza(c,1); # check norm in the wf
c0 = c # initial wavefunction for the dynamic evolution (used in the animation)

# potential for the dynamic evolution
Vpot_R = Vpot(2, z)

# saves the initial wavefunction (psi**2) and potential on a file
fv = open('initial2.dat', 'w')
format_v = "%.2f \t %.2f \t %.12g \n"
for i in range(0,Npoint-1):
    fv.write(format_v %(z[i],Vpot_R[i],np.abs(psi[i])**2))
fv.close()
Vpot_R = changeFFTposition(Vpot_R,Npoint,1)
V = changeFFTposition(Vpot_R,Npoint,0)

# evolution
t0=0.0
c = evolution(t0, Dtr, z, c0, Vpot_R, V, Ekin_K, write_evolution,plots)
