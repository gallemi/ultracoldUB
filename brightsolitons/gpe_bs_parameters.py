# coding: utf-8
import numpy as np
from gpe_bs_utilities import *
pi=np.pi


# Data block
# ------------------------------------------------------------------------------

# Scattering length

print("Interaction strenght of the soliton:")
print("\t(1) gn = -0.2\t(2) gn = -0.4\t(3) gn = -0.8")

while True:
  try:
     int_str = int(raw_input("\tWrite 1, 2 or 3: "))
  except ValueError:
     continue
  else:
     if int_str != 1 and int_str != 2 and int_str != 3:
         continue
     else:
         break

if int_str == 1:
    a_s = -0.005e-3
elif int_str == 2:
    a_s = -0.01e-3
elif int_str == 3:
    a_s = -0.02e-3

Zmax = 2.0**7                       # Grid half length
Npoint = 2**10                      # Number of grid points (min. 2**8)
Nparticle = 20e+3                   # Number of particles
whoz = 0.01                         # harmonic oscilator angular frequency
Ntime_out = 50                      # number of time steps for intermediate outputs
Dtr = 1.0e-3                        # real time step (min. 1.0e-2)
Dti = 1.0e-1                        # imaginary time step
time_final = 20.0                   # final time
Ntime_fin = int(time_final/Dtr)     # total number of time steps


# Derived quantities and parameters of the initial wavefunction (bright soliton)
# ------------------------------------------------------------------------------

# initial velocity

print("Initial velocity of the soliton:")

while True:
    try:
        velocity = int(raw_input("\tWrite an integer between 0 and 10: "))
    except ValueError:
        continue
    else:
        if velocity < 0 or velocity > 10:
            continue
        else:
            break

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
v = float(velocity)                 # initial velocity of the soliton


# Potential walls and barrier
# ------------------------------------------------------------------------------

# potential walls

print("Confine the soliton in a square box:")

while True:
    try:
        def_w = str(raw_input("\tY/N: "))
    except ValueError:
        continue
    else:
        if def_w != "Y" and def_w != "N" and def_w != "y" and def_w != "n":
            continue
        else:
            break

if def_w == "Y" or def_w == "y":
    wall = 1
elif def_w == "N" or def_w == "n":
    wall = 0

wall_h = 100.0                      # height of the walls

# external potential

print("External potential:")
print("\t(0) no potential\t(1) harmonic trap\t(2) square barrier")

while True:
    try:
        V_ext = int(raw_input("\tWrite 0, 1 or 2: "))
    except ValueError:
        continue
    else:
        if V_ext < 0 or V_ext > 2:
            continue
        else:
            break

def height_barrier(a):
    """Height of the barrier, defined as 'a' times the initial kinetic energy of the soliton
    with initial velocity 'v', defined as Ekin = 0.5 * v**2.
    Returns Ekin * a.
    """
    global v
    height = a * v**2 * 0.5
    return height

# potential barrier

if V_ext == 2:

    xb = 0.0*Zmax                       # position of the potential barrier

    print("Width of the potential barrier as function of the healing length 'xi':")
    print("\t(1) 0.5*xi\t(2) 1.0*xi\t(3) 2.0*xi")

    while True:
        try:
            barrier_w = int(raw_input("\tWrite 1, 2 or 3: "))
        except ValueError:
            continue
        else:
            if barrier_w < 1 or barrier_w > 3:
                continue
            else:
                break

    if barrier_w == 1:
        wb = 0.5 * xi
    elif barrier_w == 2:
        wb = 1.0 * xi
    elif barrier_w == 3:
        wb = 2.0 * xi

    xbr = xb + 0.5*wb                   # position of the right wall (barrier)
    xbl = xb - 0.5*wb                   # positon of the left wall (barrier)

    print("Height of the potential barrier as function of the kinetic energy 'Ec':")

    while True:
        try:
            barrier_h = float(raw_input("\tWrite a float between 0 and 2: "))
        except ValueError:
            continue
        else:
            if barrier_h < 0.0 or barrier_h > 2.0:
                continue
            else:
                break

    hb = height_barrier(barrier_h)

else:
    xb, hb, wb = 0, 0, 0
