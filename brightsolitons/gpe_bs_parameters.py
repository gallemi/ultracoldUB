# coding: utf-8
import numpy as np
from gpe_bs_utilities import *
pi=np.pi

# Data block
# ------------------------------------------------------------------------------

Zmax = 2.0**7                       # Grid half length
Npoint = 2**10                      # Number of grid points (min. 2**8)
Nparticle = 20e+3                   # Number of particles
a_s = -0.01e-3                      # scattering length (cases: -0.005, -0.01, -0.02)
whoz = 0.01                         # harmonic oscilator angular frequency
Ntime_out = 20                      # number of time steps for intermediate outputs
Dtr = 1.0e-3                        # real time step (min. 1.0e-2)
Dti = 1.0e-1                        # imaginary time step
time_final = 20.0                   # final time
Ntime_fin = int(time_final/Dtr)     # total number of time steps


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
x0 = -50.0                          # initial position of the soliton
v = +6.0                            # initial velocity of the soliton (from 0 to 10)


# Potential walls and barrier
# ------------------------------------------------------------------------------

def height_barrier(a):
    """Height of the barrier, defined as 'a' times the initial kinetic energy of the soliton
    with initial velocity 'v', defined as Ekin = 0.5 * v**2.
    Returns Ekin * a.
    """
    global v
    height = a * v**2 * 0.5
    return height

xb = 0.0*Zmax                       # position of the potential barrier
hb = height_barrier(0.5)            # height of the potential barrier
wb = 0.5 * xi                       # width of the potential barrier
wall_h = 100.0                      # height of the walls, fixed
wall = 0                            # defines the walls of the box if not 0
xbr = xb + 0.5*wb                   # position of the right wall (barrier)
xbl = xb - 0.5*wb                   # positon of the left wall (barrier)
