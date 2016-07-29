
# definitivo (o eso creo)
# coding: utf-8

# ## FFT solver for 1D Gross-Pitaevski equation

# We look for the complex function $\psi(x)$ satisfying the GP equation
#
# $ i\partial_t \psi = -\frac{1}{2}(i\partial_x - \Omega)^2\psi+ V(x)\psi + g|\psi|^2\psi $,
#
# with periodic boundary conditions.
#
# Integration: pseudospectral method with split (time evolution) operator;
# that is evolving in real (R) or momentum (K) space according to the operators
# in the Hamiltonian, i.e.
# first we evaluate
#
# $\hat{\psi}(x,\frac{t}{2})=\cal{F}^{-1}\left[\exp\left(-i \frac{\hbar^2 k^2}{2} \frac{t}{2}\right)\,\psi(k,0)\right] $
#
# and later
#
# $\psi(k,t) = \exp(-i \frac{\hbar^2 k^2}{2} \frac{t}{2})\,
# \cal{F}\left[\exp\left(-i (V(x)+|\hat{\psi}(x,\frac{t}{2})|^2)\, t \right)\,\hat{\psi}(x,\frac{t}{2}) \,
# \right]$
#
# where $\cal{F}$ is the Fourier transform.
#
# Program: Evolution in real time. We evolve in real time the wave pack, the user decides the type of movement (free or osc. harm.).
# The program uses local folders called "gpe_fft_utilities" and "wave_functions.py". There is different type of plots and an important function called
# "changeFFTposition" that it is important to Fourier transform.
# _______________________________________________________________________________
# _______________________________________________________________________________

# Import libraries and general definitions:
#_________________________________________________________________________________________________

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
from gpe_fft_utilities import * # local folder utilities
from wave_functions import *
import numpy.linalg as lin
from pylab import*
import os, glob

close('all')
pi=np.pi


# Data block
# ________________________________________________________________________________________________

# In[2]:

# User decides the time
while True:
    try:
        osci =float(raw_input('introduce el numero de unidades de tiempo que dura la simulacion entre 1 y 10'))
        while ((osci<1) or (osci>10)):
            print "ERROR: la simulacion debe de estar entre un rango de 1 y 10"
            osci=float(raw_input("introduce el numero de unidades de tiempo que dure la simulacion"))
        break
    except ValueError:
        print("Escribe un numero")

Zmax = 50.0              # Grid half length
Npoint =512              # Number of grid points
Nparticle = 500          # Number of particles
a_s = 0.0                # scattering length

# User decides the type of movement (free movement or harm. osc.)
while True:
    try:
        harm=float(raw_input('introduce si quieres movimiento oscilatorio harmonico o no (1=si ; 0=no)'))
        while (harm != 1 and harm!=0):
            print 'ERROR: introduce 0 o 1'
            harm=float(raw_input('introduce 1 o 0'))
        break
    except ValueError:
        print ("Escoge las opciones que se te han dado")

if (harm==1):
    whoz = 1.0               # harmonic oscilator angular frequency
if (harm==0):
    whoz = 0.0               # harmonic oscilator angular frequency

Omega = pi/(2*Zmax)          # reference frame velocity
Dtr = 1.0e-3                 # real time step
Dti = 1.0e-3                 # imaginary time step
Ntime_fin=int(osci*1000)     # total number of time steps
Ntime_out = 100              # number of time steps for intermediate outputs

# We choose the initial position of wave pack:

while True:
    try:
        x0=float(raw_input("introduce posicion inicial del paquete de 0 a 5"))
        while (np.abs(x0)>(5)): # tolerance for the initial position of the soliton
            print "ERROR: la posicion inicial del paquete debe estar dentro del rango de 0 a 5"
            x0=float(raw_input("introduce posicion inicial del soliton"))
        break
    except ValueError:
        print ("Escoge un numero")

# Print evolution data:

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
# _________________________________________________________________________________________

# In[3]:

NormWF = 1.0/(2*Zmax)           # Wave function (WF) norm
gint = 2*a_s*Nparticle*NormWF   # interaction (nonlinear-term) strength
Dz = 2*Zmax/Npoint              # length step size
Dk = pi/Zmax                    # momentum step size
Kmax = Dk*(Npoint//2)           # maximum momentum
Dt = Dtr-1j*Dti                 # complex time
Ninter = Ntime_fin/Ntime_out    # Number of outputs with the intermediate states
print(" Characteristic interaction energy = %g"%(gint))


# Grid definitions: physical and momentum space
# __________________________________________________________________________________________

# In[5]:

z = np.arange(-Zmax+Dz,Zmax+Dz,Dz)  # physical (R-space) grid points in ascending order
zp = changeFFTposition(z,Npoint,1)  # (R-space) grid points with FFT order

kp = np.arange(-Kmax+Dk,Kmax+Dk,Dk) # physical (K-space) grid points in ascending order
kp = changeFFTposition(kp,Npoint,1) # (K-space) gridd points with FFT order


# Define kinetic energy and potential:
# ________________________________________________________________________________________

# In[6]:

Ekin_K = 0.5*(kp**2) # Kinetic energy in K space

# Potential energy in R space:
# Harmonic oscillator with angular frequency whoz:

Vpot_R = 0.5*whoz*zp**2


# Main functions:
# ________________________________________________________________________________________

# In[7]:

def Energy(c): # Energy (per particle) calculation
    global gint, Vpot_R, Ekin_K, Npoint
    ek = sum(Ekin_K*abs(c)**2)              # Kinetic energy in K
    psi = ifft(c)*Npoint;                   # wf FFT to R
    ep = sum(Vpot_R*abs(psi)**2)/Npoint;    # Potential energy
    ei = 0.5*gint*sum(abs(psi)**4)/Npoint;  # Interaction energy
    em =  ek+ep+ei;                         # average energy
    chem_pot = em+ei;                       # chemical potential
    return em, chem_pot, ek, ep, ei

def T_K (t,Dt,psi): # Action of the time evolution operator over state c in K space
    global Ekin_K
    #psi is the wave function in K space
    # t is the time (which is not used for time independant Hamiltonians)
    # Dt is the complex time step

    return np.exp(-1j*0.5*Dt*Ekin_K)*psi # return action on psi

def T_R_psi(t,Dt,psi): # Action of the time evolution operator over state c in R space
    global gint, Vpot_R
    # Includes the external potential and the interaction operators:
    #       T_R_psi = exp(-i Dt (Vpot_R+ gint|psi|^2) ) c
    # psi is the wave function in R space
    # t is the time (which is not used for time independant Hamiltonians)
    # Dt is the complex time step

    return np.exp( -1j*Dt*(Vpot_R + gint*(abs(psi)**2)) )*psi # return action on psi

def normaliza(c): # normalization to 1
    norm = lin.norm(c)
    if ((norm-1.0)>1.0e-4): # check norm
        print("normalization from: ",norm)
    return c/norm


# Plots of the initial state:
#____________________________________________________________________________________________

# In[9]:

t=0.0
c=normaliza(gaussian(zp,Npoint,x0,0,1.0,0.0))
cc = ifft(c)*Npoint*NormWF**0.5      # FFT from K3 to R3 and include the wf norm
psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function

#plot different propieties of psi:

plot_density(z,psi,Zmax,t)
plot_phase(z,psi,Zmax,t)
plot_real_imag(z,psi,Zmax,t)

print("Energies:          Emed    mu    Ekin    Epot    Eint")
print("         initial = %g %g %g %g %g"%(Energy(c)))

psi_sol=psi                  # we chance name variable

# Choose initial wave function and evolve in real time:
# __________________________________________________________________________________________

# In[11]:
# temporal evolution soliton

psi_sol=changeFFTposition(psi_sol,Npoint,1)
c0=normaliza(fft(psi_sol)/Npoint) # initial wave function

print("Energies in evolution real time:          Emed    mu    Ekin    Epot    Eint")
print("         initial = %g %g %g %g %g"%(Energy(c0)))

# evolution in time: parameters
t0=0.0
c=c0
j=0
t=0
tevol=np.empty([Ninter+1])     # time vector
pos_minus=np.empty([Ninter+1]) # put the minus position in a vector
val_minus=np.empty([Ninter+1]) # put the minus value in a vector
energi=np.empty([5])           # put the energies in a vector

tevol[0]=t0

# where the files of evolution will be saved
dir = "wp_evolution" # name of the directory
if (not os.path.exists(dir)): # creates the directory
    os.makedirs(dir)
else: # removes its contents if it already exists
    print "Directory %r already exists. Do you want to continue?" % (dir)
    while True:
        try:
            ans = str(raw_input("\tY/N: ")) # pause, answer to continue
        except ValueError:
            continue
        else:
            if (ans!='Y' and ans!='y' and ans!='N' and ans!='n'):
                continue
            else:
                break
    if (ans == 'N' or ans == 'n'):
        sys.exit("End of program")


# Open files
file=open('energies.txt','w')
file.write('Tabla donde se muestran diversos valores de la energia a lo largo del movimiento del soliton.\n')
file.write('Tiempo\tEnergia media\tPotencial quimico\tEnergia cinetica\tEnergia potencial\tEnergia interna\n' )

f4=plt.figure()
for i in range(1, Ntime_fin+1): # time evolution cicle
    t += Dt.real
    psi=ifft(T_K(t0,Dt.real,c))*Npoint
    c=T_K(t0,Dt.real,fft(T_R_psi(t0,Dt.real,psi))/Npoint)
    c = normaliza(c); # check norm in the wf


    if(not(i%Ntime_out)):
        j+=1
        tevol[j] = t
# Write energies from function Energy
        energi=(Energy(c))
        file.write('%s\t' %t)
        file.write('%g\t%g\t%g\t%g\t%g\n' %(Energy(c)))
# Representation of intermediate solutions
        cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
        psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function

        plt.title('Evolution in time'%(tevol[Ninter]),fontsize=15)
        plt.xlabel('$x/a_{ho}$',fontsize=15)
        plt.xticks(np.arange(-Zmax, Zmax+1,Zmax/2))
        plt.locator_params('y',nbins=3)
        plt.plot(z, abs(psi)**2, 'b--',label='$|\psi|^2$') # plot density
#        plt.plot(z, psi.real, 'r.',label='real$(\psi)$')
#        plt.plot(z, psi.imag, 'b--',label='imag$(\psi)$')
#        plt.plot(z, np.angle(psi), 'b.',label='$Arg(\psi)$')
        f4.show()
        psi*=np.exp(1j*pi/3) # This is useful to plot the wave function phase.

# Writes wave function
        file2=open('./%s/WfWd-%08d.txt'%(dir, j),'w')
        file2.write('Tiempo=%s\n' %(t))
        file2.write('Datos de interes: N.particulas=%g\tPar.Interaccion=%g\tLong.caja=%g\tN.puntos=%g\tFreq.Oscilador=%g\tPot. quim.=%s\n ' %(Nparticle,gint,2*Zmax,Npoint,whoz,energi[1]))
        file2.write('x\tDensidad\tFase\tRe\tIm\tV(x)\n')
        for i in range (0,int(2*Zmax/Dz)):
            file2.write("%s\t%s\t%s\t%s\t%s\t%s \n" %(z[i],(abs(psi)**2)[i],(np.angle(psi))[i],psi.real[i],psi.imag[i],changeFFTposition(abs(c)**2,Npoint,0)[i]))

plt.show()
file.close()
file2.close()
# Prints final energy
print("         final = %g %g %g %g %g"%(Energy(c)))
