# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as lin
from scipy.fftpack import fft, ifft
import scipy.optimize as opt
from sys import exit
import os, glob
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
       wb      halfwidth of the barrier
       hb      height of the barrier
       wall    height of the walls
    """
    global Zmax, Npoint, whoz, xb, wb, hb, wall
    Vpot_R = np.empty([Npoint]) # vector of the same length as z, initially empty
    m = 1000                    # slope used to define the walls and barrier
    wall_x = 0.5*Zmax           # position of the walls relative to the halfwidth of the box
    # defines the walls of the box
    if wall != 0.0 :
        Vpot_R = wall - wall/(1.0+np.exp(m*(z-wall_x))) + wall/(1.0+np.exp(m*(z+wall_x)))
    else:
        Vpot_R = 0.0
    # defines the potential (barriers, harmonic trap, etc.)
    if(n==0):
        Vpot_R = Vpot_R
    elif(n==1):
        Vpot_R = 0.5*whoz**2*z**2
    elif(n==2):
        Vpot_R = Vpot_R + hb/(1.0+np.exp(m*(z-xb-wb))) - hb/(1.0+np.exp(m*(z-xb+wb)))
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

# Evolution
# ------------------------------------------------------------------------------

def evolution(t0, Dt, z, c0, Vpot_R, V, Ekin_K, write_ev):
    """Calculates the evolution of the wavefunction c.
       t0         initial time
       Dt         time step (either imaginary or real)
       z          R-space grid points
       c0         initial wavefunction (FFT order)
       Vpot_R     external potential (FFT order)
       V          external potential (physical order)
       Ekin_K     kinetic energy (FFT order)
       write_ev   writes data into files if 0
    Global variables:
       Ntime_out  number of time steps for intermediate outputs
       Ntime_fin  total number of time steps
       gint       interaction strength * NormWF
       Npoint     number of points in the grid
       NormWF     norm of the wavefunction
       Zmax       halfwidth of the grid
       Nparticle  number of particles
       a_s        scattering length
       x0         initial position of the soliton
       v          initial velocity of the soliton
       whoz       harmonic oscillator frequency
       xb         position of the barrier
       wb         halfwidth of the barrier
       hb         height of the barrier
       wall       height of the walls
    Returns the final wavefunction.
    If write_ev == 0:
       the wavefunction (psi) and related magnitudes are written on a file for
           a certain time steps (a different file for each step).
       the file name includes the corresponding time x 1000.
       such files are created on a directory inside the current directory.
       the program checks if the directory where the files will be stored exists
           and if so deletes its contents (otherwise creates a new directory).
    """
    global Ntime_out, Ntime_fin, gint, Npoint, NormWF, Zmax, Nparticle, a_s, x0, v, whoz, xb, wb, hb, wall
    # where the files of evolution will be saved
    if (write_ev==0):
        dir = "bright_barrier-v%s-x0%s-xb%s" % (v, x0, xb) # name of directory
        name = "WfBs" # general name of the files
        if (not os.path.exists(dir)): # creates the directory
            os.makedirs(dir)
        else: # removes its contents if it already exists
            print "Directory %r already exists..." % (dir)
            raw_input(" ") # pause, any input will work
            files = "%s/%s_*.dat" % (dir, name)
            for f in glob.glob(files):
                os.remove( f )

    # open files for normalization, energies, and some values of c
    if isinstance(Dt, complex):
        print "(Imaginary time evolution)"
        fn = open('normalization_imag.dat', 'w')
        fe = open('energies_imag.dat', 'w')
        fc = open('c_imag.dat', 'w')
    else:
        print "(Real time evolution)"
        fn = open('normalization_real.dat', 'w')
        if (write_ev==0):
            fe = open('./%s/energies.dat' % dir, 'w')
        else:
            fe = open('energies_real.dat', 'w')
        fc = open('c_real.dat', 'w')

    # wavefunction and counters
    Ninter = Ntime_fin//Ntime_out # number of outputs (intermediate states)
    c = c0
    j=0; t=t0; l=0

    # define vectors for the energies and time
    tevol=np.empty([Ninter+1])
    energy_cicle=np.empty([Ninter+1,5])

    # headers and formats for files
    format_c = "%.2f \t %.12g \t %.12g \t %.12g \n"
    format_e = "%.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \n"
    fe.write("# %s\t%s\t%s\t%s\t%s\t%s\n" %("time","total energy","chemical potential","kinetic energy","potential energy","interaction energy"))
    header_variables = "x", "|psi|^2", "phase", "Re(psi)", "Im(psi)", "V(x)", Nparticle, a_s, 2*Zmax, Npoint, Ntime_fin, Dt, x0, v, xb, 2*wb, hb, wall, energy_cicle[j,1], energy_cicle[j,0]
    header_format = "# %s"+("\t%s")*5+"\n" + "# Number of particle = %s; scattering length = %s; Grid width = %s; Number of points = %s; Total number of time steps = %s; time step = %s " + "initial position of the soliton = %s; initial velocity of the soliton = %s; position of the barrier = %s; width of the barrier = %s; height of the barrier = %s; height of the potential walls = %s" + " results: chemical potential = %s; total energy = %s \n"
    format_psi = "%.2f" + ("\t %.12g")*5 + "\n"

    # energy and time at t0 (prints and writes)
    tevol[0]=t0
    energy_cicle[0,:] = Energy(c, Vpot_R, Ekin_K)
    fe.write(format_e %(tevol[j], energy_cicle[j,0], energy_cicle[j,1], energy_cicle[j,2], energy_cicle[j,3], energy_cicle[j,4]))
    print("Energies:          Emed    mu    Ekin    Epot    Eint")
    print("         initial = %g %g %g %g %g"%(Energy(c0, Vpot_R, Ekin_K)))

    # wavefunction
    cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
    psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function

    # plots initial state and prepares plot of intermediate states
    f4=plt.figure()
    if isinstance(Dt, complex):
        plt.title('Evolution of the initial wavefunction (imaginary time)',fontsize=15)
    else:
        plt.title('Evolution of the initial wavefunction (real time)',fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    plt.xticks(np.arange(-Zmax, Zmax+1,Zmax/2))
    plt.axis((-Zmax,Zmax,0,0.3))
    plt.locator_params('y',nbins=3)
    plt.plot(z, abs(psi)**2, 'r-',label='$|\psi_0|^2$')
    plt.plot(z, changeFFTposition(Vpot_R,Npoint,0), 'g--',label='potential')
    plt.legend(fontsize=15)

    # writes initial state on the corresponding file
    if (write_ev==0):
        fpsi = open('./%s/%s-%08d.dat' %(dir,name,tevol[j]*1000), 'w')
        fpsi.write(header_format %(header_variables))
        for k in range(0,Npoint-1):
            fpsi.write(format_psi %(z[k], np.abs(psi[k]**2), np.angle(psi[k]), psi[k].real, psi[k].imag, V[k]))

    # time evolution cicle
    for i in range(1, Ntime_fin+1):
        t += np.abs(Dt)
        psi=ifft(T_K(Dt, Ekin_K)*c)*Npoint
        c=T_K(Dt, Ekin_K)*fft(T_R_psi(t0,Dt,psi,Vpot_R))/Npoint
        fc.write(format_c % (t , np.abs(c[int(Npoint/2-0.1*Npoint)]) , np.abs(c[int(Npoint/2)]) , np.abs(c[int(Npoint/2+0.1*Npoint)]))) # write some values of c on a file
        c = normaliza(c,fn); # check norm in the wf
        V = changeFFTposition(Vpot_R,Npoint,0) # in physical order
        if(not(i%Ntime_out)):
            j+=1
            tevol[j] = t
            energy_cicle[j,:] = Energy(c, Vpot_R, Ekin_K)
            fe.write(format_e %(tevol[j], energy_cicle[j,0], energy_cicle[j,1], energy_cicle[j,2], energy_cicle[j,3], energy_cicle[j,4]))
            if(not(i%100)):
                cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
                psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function
                if(not(i%500)):
                    plt.plot(z, abs(psi)**2, 'b--',label='$|\psi|^2$') # plot density
                    if(l==0):
                        plt.legend(fontsize=15)
                    l+=1
                if(write_ev==0):
                    fpsi = open('./%s/%s-%08d.dat' %(dir,name,tevol[j]*1000), 'w')
                    fpsi.write(header_format %(header_variables))
                    format_psi = "%.2f" + ("\t %.12g")*5 + "\n"
                    for k in range(0,Npoint-1):
                        fpsi.write(format_psi %(z[k], np.abs(psi[k]**2), np.angle(psi[k]), psi[k].real, psi[k].imag, V[k]))
                    fpsi.close()

    # prints final energies
    print("         final   = %g %g %g %g %g"%(Energy(c, Vpot_R, Ekin_K)))
    print("Energy change at last step  = %g"%(energy_cicle[Ninter,0]-energy_cicle[Ninter-1,0]))
    f4.show(); plt.show()

    # plots energies
    plot_convergence(tevol,energy_cicle[:,0],energy_cicle[:,1],energy_cicle[:,2],energy_cicle[:,3],energy_cicle[:,4],Ninter); plt.show()

    # other plots (for the final state)
    plot_phase(z,psi,Zmax,t); plt.show()
    plot_real_imag(z,psi,Zmax,t); plt.show()

    # closes files
    fn.close(); fe.close(); fc.close()

    return c


# Other utilities and plots
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

def plot_convergence(x,y0,y1,y2,y3,y4,n):
    f1=plt.figure()
    plt.title('Convergence',fontsize=15)
    plt.xlabel('time ($t \, \\omega_{ho}$)',fontsize=15)
    plt.ylabel('Energy per particle ($E/\\hbar \,\\omega_{ho}$)',fontsize=15)
    plt.xticks(np.arange(0, x[n]+1,x[n]/5))
    plt.locator_params('y',nbins=3)
    plt.plot(x, y0, 'r-',label="$E_{med}$")
    plt.plot(x, y2, 'b-',label="$E_{kin}$")
    plt.plot(x, y3, 'g-',label="$E_{pot}$")
    plt.plot(x, y4, 'y-',label="$E_{int}$")
    plt.plot(x, y1, 'm',label="$\\mu$")
    plt.legend(fontsize=15)
    f1.show()

def plot_density(z,psi,Lz,t):
    f2=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, abs(psi)**2, 'b--',label='$|\psi|^2$') # plot density
    plt.legend(fontsize=15)
    f2.show()

def plot_phase(z,psi,Lz,t):
    f3=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, np.angle(psi), 'b.',label='$Arg(\psi)$')
    plt.legend(fontsize=15)
    f3.show()

def plot_real_imag(z,psi,Lz,t):
    f3=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, psi.real, 'r.',label='real$(\psi)$')
    plt.plot(z, psi.imag, 'b--',label='imag$(\psi)$')
    plt.legend(fontsize=15)
    f3.show()
