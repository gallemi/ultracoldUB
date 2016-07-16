# coding: utf-8
import numpy as np
from gpe_bs_utilities import *
pi=np.pi

# Data block
# ------------------------------------------------------------------------------

Zmax = 2.0**6               # Grid half length
Npoint = 2**10              # Number of grid points
Nparticle = 20              # Number of particles
a_s = -0.01                 # scattering length
whoz = 0.01                 # harmonic oscilator angular frequency
Omega = 0 * pi/(2*Zmax)     # reference frame velocity
Ntime_out = 20              # number of time steps for intermediate outputs
Dtr = 1.0e-3                # real time step
Dti = 1.0e-2                # imaginary time step
Ntime_fin = 10000           # total number of time steps


# Derived quantities and parameters of the initial wavefunction (bright soliton)
# ------------------------------------------------------------------------------

NormWF = 1.0/(2*Zmax)               # Wave function (WF) norm
Dz = 2*Zmax/Npoint                  # length step size
Dk = pi/Zmax                        # momentum step size
Kmax = Dk*(Npoint//2)               # maximum momentum
Dt = Dtr-1j*Dti                     # complex time
Ninter = Ntime_fin//Ntime_out       # Number of outputs with the intermediate states
gn = 2*a_s*Nparticle                # interaction (nonlinear-term) strength
gint = gn*NormWF                    # int. strenght, with factor NormWF
xi = 1.0/(np.abs(gn)**2*0.5)**0.5   # healing length
x0 = 0.0                            # initial position of the soliton
wall = 100.0                        # height of the walls of the box
xb = 0.3*Zmax                       # position of the potential barrier
hb = 4.5                            # height of the potential barrier
wb = 0.0575*xi                      # halfwidth of the potential barrier
v = +3.0 * 2.0*Zmax                 # initial velocity of the soliton
