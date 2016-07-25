# coding: utf-8
import numpy as np
import numpy.linalg as lin
from scipy.fftpack import fft, ifft
from sys import exit
from gpe_bs_parameters import *


# Initial wavefunction (bright soliton) and potential
# ------------------------------------------------------------------------------

def bright_soliton(x,n,x0,gn):
    """Bright soliton: 1/cosh(x/xi) wf, analytic solution
       x  = grid points
       n  = number of grid points
       x0 = soliton position
       gn = (negative) interaction energy or double of the chemical potential
    """
    if (gn >= 0):
        exit(" error: interaction must be negative...")
    psi0 = 0.5*np.sqrt(np.abs(gn))
    xi = 1.0/(2.0*np.abs(gn)*np.abs(psi0)**2)**0.5
    fx = psi0 * 1.0/np.cosh((x-x0)/(xi * np.sqrt(2.0)))   # wave function in R3
    return fft(fx)/n    # FFT to K3

def Vpot(n, z):
    """External potential according to n. Cases:
       n=0     no potential
       n=1     harmonic trap
       n=2     square potential barrier
    The walls of the box containing the soliton are defined allways unless
    the height of the walls is set to zero.
    Parameters (introduced as arguments):
       z       R-space grid points
    Global variables:
       Zmax    halfwidth of the grid
       Npoint  number of grid points
       whoz    harmonic oscillator frequency
       xb      position of the barrier
       xbr     position of the right wall (barrier)
       xbl     position of the left wall (barrier)
       hb      height of the barrier
       wall    defines potential walls if not 0
    """
    global Zmax, Npoint, whoz, xb, xbr, xbl, hb, wall, wall_h
    Vpot_R = np.empty([Npoint]) # vector of the same length as z, initially empty
    m = 1000                    # slope used to define the walls and barrier
    wall_x = 0.75*Zmax          # position of the walls relative to the halfwidth of the box
    # defines the walls of the box
    if wall != 0 :
        Vpot_R = wall_h - wall_h/(1.0+np.exp(m*(z-wall_x))) + wall_h/(1.0+np.exp(m*(z+wall_x)))
    else:
        Vpot_R = 0.0*z
    # defines the potential (barriers, harmonic trap, etc.)
    if(n==0):
        Vpot_R = Vpot_R
    elif(n==1):
        Vpot_R = 0.5*whoz**2*z**2
    elif(n==2):
        Vpot_R = Vpot_R + hb/(1.0+np.exp(m*(z-xbr))) - hb/(1.0+np.exp(m*(z-xbl)))
    else:
        Vpot_R = Vpot_R
    return Vpot_R


# operators and energies
# ------------------------------------------------------------------------------

def Energy(c, Vpot_R, Ekin_K):
    """Energy (per particle) calculation for a certain time.
       c         initial wavefunction (FFT ordered)
       Vpot_R    potential as function of the position
       Ekin_K    kinetic energy as function of the position
    Global variables:
       gint      interaction strength * NormWF (norm)
       Npoint    number of points in the grid
    Energy terms:
       ek        kinetic energy
       ep        potential energy
       ei        interaction energy
       em        average energy
       chem_pot  chemical potential
    Returns the energy terms.
    """
    global gint, Npoint
    ek = sum(Ekin_K*abs(c)**2)              # Kinetic energy in K
    psi = ifft(c)*Npoint;                   # wf FFT to R
    ep = sum(Vpot_R*abs(psi)**2)/Npoint;
    ei = 0.5*gint*sum(abs(psi)**4)/Npoint;
    em =  ek+ep+ei;
    chem_pot = em+ei;
    return em, chem_pot, ek, ep, ei

def T_K(Dt, Ekin_K):
    """Time evolution operator in K space (for second order accuracy).
       Dt      complex time step
       Ekin_K  kinetic energy (function of the position)
    """
    T_K = np.exp(-1j*0.5*Dt*Ekin_K)
    return T_K

def T_R_psi(t, Dt, psi, Vpot_R):
    """Action of the time evolution operator over state c in R space
       t       time (not used for a time independent Hamiltonian)
       Dt      complex time step
       psi     wavefunction in R space
       Vpot_R  external potential
    Global variables:
       gint    interaction strength * NormWF (norm)
    Includes the external potential and the interaction operators:
       T_R_psi = exp(-i Dt (Vpot_R+ gint|psi|^2) ) c
    Returns the action on psi.
    """
    global gint
    return np.exp( -1j*Dt*(Vpot_R + gint*(abs(psi)**2)) )*psi



# Other utilities
# ------------------------------------------------------------------------------

def normaliza(c,file):
    """Checks normalization of the wavefunction c and sets norm to 1.
    Norm is written on a file unless variable 'file' is either 0 or 1."""
    frmt = "%.10f \n"
    norm = lin.norm(c)
    if file != 0 and file != 1: #and != 1:
        file.write(frmt %(norm))
    elif file == 1:
        print("Initial norm: %g" %(norm**2) )
    return c/norm

def changeFFTposition(f,N,j):
    """Change the order in vectors from FFT
       f(0...N-1) is the vector to order
       N is the vector dimension
       j is a switch indicating the change direction
       physical order is f = [(-(Zmax-Dz):Dz:-Dz) (0:Dz:Zmax)]
       FFT order is f = [(0:Dk:kmax) (-(kmax-kz):kz:-kz)]
    """
    f1 = f*1
    if (j==1):  # from physical to FFT order
        for i in range(0,N//2-1) :
            f1[i] = f[N//2-1+i];
            f1[N//2+1+i] = f[i];
        f1[N//2-1] = f[N-2];
        f1[N//2] = f[N-1];
    elif (j==0): # from FFT to physical order
        for i in range(0,N//2-1) :
            f1[i] = f[N//2+1+i];
            f1[N//2+1+i] = f[i+2];
        f1[N//2-1] = f[0];
        f1[N//2] = f[1];
    else:
        print("error in changeFFTposition(f,N,j): j must be 0 or 1...")
    return f1

def integral(x,dz,z,i,f):
    """Integrates function "x" from i to f.
    x is a function (vector)
    dz is the length step size
    z is the grid
    i is the initial piont
    f is the last point
    k, j are counters
    """
    k = 0
    if len(x) != len(z):
        print "Error: the function has not the same length as the grid"
    inte = 0.0
    for j in z:
        if j>=i and j<=f:
            inte += x[k]*dz
        k += 1
    return inte

def list_integrals(x,z):
    global Dz, Zmax, xbl, xbr
    ileft = integral(x,Dz,z,-Zmax,xbl)
    iin = integral(x,Dz,z,xbl,xbr)
    iright = integral(x,Dz,z,xbr,Zmax)
    return ileft, iin, iright
