# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft, ifft
import os, glob
from gpe_bs_parameters import *
from gpe_bs_plots import *
from gpe_bs_utilities import *
import sys


# Evolution
# ------------------------------------------------------------------------------

def evolution(t0, Dt, z, c0, Vpot_R, V, Ekin_K, write_ev, plots):
    """Calculates the evolution of the wavefunction c.
       t0         initial time
       Dt         time step (either imaginary or real)
       z          R-space grid points
       c0         initial wavefunction (FFT order)
       Vpot_R     external potential (FFT order)
       V          external potential (physical order)
       Ekin_K     kinetic energy (FFT order)
       write_ev   writes data into files if 0
       plots      plots data if 0
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
       wb         width of the barrier
       hb         height of the barrier
       wall_h       height of the walls
    Returns the final wavefunction.
    If write_ev == 0:
       the wavefunction (psi) and related magnitudes are written on a file for
           a certain time steps (a different file for each step).
       the file name includes the corresponding time x 1000.
       such files are created on a directory inside the current directory.
       the program checks if the directory where the files will be stored exists
           and if so deletes its contents (otherwise creates a new directory).
    """
    global Ntime_out, Ntime_fin, gint, Npoint, NormWF, Zmax, Nparticle, a_s, x0, v, whoz, xb, wb, hb, wall_h, wall
    # where the files of evolution will be saved
    if (write_ev==0):
        dir = "bs-gn_%0d-v_%0d-hb_%0d-wb_%0d" % (np.abs(gn)*10, np.abs(v), np.abs(hb), np.abs(wb)) # name of directory
        name = "WfBs" # general name of the files
        if (not os.path.exists(dir)): # creates the directory
            os.makedirs(dir)
        else: # removes its contents if it already exists
            print "Directory %r already exists." % (dir)
            while True:
                try:
                    ans = str(raw_input("\t Do you want to continue? (Y/N): ")) # pause, answer to continue
                except ValueError:
                    continue
                else:
                    if (ans!='Y' and ans!='y' and ans!='N' and ans!='n'):
                        continue
                    else:
                        break
            if (ans == 'N' or ans == 'n'):
                sys.exit("End of program")            
                
            
            files = "%s/%s_*.dat" % (dir, name)
            for f in glob.glob(files):
                os.remove( f )

    # open files for normalization, energies, and some values of c
    if isinstance(Dt, complex):
        print "(Imaginary time evolution)"
        fn = open('normalization_imag.dat', 'w')
        fe = open('energies_imag.dat', 'w')
        #fc = open('c_imag.dat', 'w')
        fpsi_all = open('evolution_imag.dat', 'w')
    else:
        print "(Real time evolution)"
        fn = open('normalization_real.dat', 'w')
        if (write_ev==0):
            fe = open('./%s/energies.dat' % dir, 'w')
        else:
            fe = open('energies_real.dat', 'w')
        #fc = open('c_real.dat', 'w')
        fpsi_all = open('evolution_real.dat', 'w')

    # wavefunction and counters
    Ninter = Ntime_fin//Ntime_out # number of outputs (intermediate states)
    c = c0
    j=0; t=t0; l=0

    # define vectors for the energies, time and % of the wf (matrix)
    tevol=np.empty([Ninter+1])
    energy_cicle=np.empty([Ninter+1,5])
    wave_function = np.empty([Ninter+1,3])

    # headers and formats for files
    format_e = "%.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \t %.10f \n"
    fe.write("# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("time","total energy","chemical potential","kinetic energy","potential energy","interaction energy", "left integral", "inside integral", "right integral"))
    header_variables = "x", "|psi|^2", "phase", "Re(psi)", "Im(psi)", "V(x)", Nparticle, a_s, 2*Zmax, Npoint, Ntime_fin, Dt, x0, v, xb, wb, hb, energy_cicle[j,1], energy_cicle[j,0]
    header_format = "# %s"+("\t%s")*5+"\n" + "# Number of particle = %s; scattering length = %s; Grid width = %s; Number of points = %s; Total number of time steps = %s; time step = %s " + "initial position of the soliton = %s; initial velocity of the soliton = %s; position of the barrier = %s; width of the barrier = %s; height of the barrier = %s;" + " results: chemical potential = %s; total energy = %s \n"
    format_psi = "%.2f" + ("\t %.12g")*5 + "\n"

    # energy and time at t0 (prints)
    tevol[0]=t0
    energy_cicle[0,:] = Energy(c, Vpot_R, Ekin_K)
    print("Energies:          Emed    mu    Ekin    Epot    Eint")
    print("         initial = %g %g %g %g %g"%(Energy(c0, Vpot_R, Ekin_K)))

    # wavefunction
    cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
    psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function
    wave_function[0,:] = list_integrals(np.abs(psi)**2,z)

    # writes initial energies and wavefunction on a file
    fe.write(format_e %(tevol[0], energy_cicle[0,0], energy_cicle[0,1], energy_cicle[0,2], energy_cicle[0,3], energy_cicle[0,4], wave_function[0,0], wave_function[0,1], wave_function[0,2]))

    # state at t0
    for k in range(0,Npoint-1):
        fpsi_all.write(format_psi %(z[k], np.abs(psi[k]**2), np.angle(psi[k]), psi[k].real, psi[k].imag, V[k]))

    # plots initial state and prepares plot of intermediate states
    if(plots==0):
        f4=plt.figure()
        if isinstance(Dt, complex):
            plt.title('Evolution of the initial wavefunction (imaginary time)',fontsize=15)
        else:
            plt.title('Evolution of the initial wavefunction (real time)',fontsize=15)
        plt.xlabel('$z$',fontsize=15)
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
        c = normaliza(c,fn); # check norm in the wf
        V = changeFFTposition(Vpot_R,Npoint,0) # in physical order

        if(not(i%Ntime_out)):
            j+=1
            tevol[j] = t
            cc = ifft(c)*Npoint*NormWF**0.5 # FFT from K3 to R3 and include the wf norm
            psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function
            energy_cicle[j,:] = Energy(c, Vpot_R, Ekin_K)
            wave_function[j,:] = list_integrals(np.abs(psi)**2,z)
            fe.write(format_e %(tevol[j], energy_cicle[j,0], energy_cicle[j,1], energy_cicle[j,2], energy_cicle[j,3], energy_cicle[j,4], wave_function[j,0], wave_function[j,1], wave_function[j,2]))

            # plots and writes
            if(not(i%100)):
                # file for the animation with gnuplot
                for k in range(0,Npoint-1):
                    fpsi_all.write(format_psi %(z[k], np.abs(psi[k]**2), np.angle(psi[k]), psi[k].real, psi[k].imag, V[k]))
                fpsi_all.write("\n \n")

                # plots density
                if(plots==0):
                    if(not(i%500)):
                        plt.plot(z, abs(psi)**2, 'b--',label='$|\psi|^2$')
                        if(l==0):
                            plt.legend(fontsize=15)
                        l+=1

                # writes a file for each timestep
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
    print("  E(final) - E(initial) = %g"%(np.abs(energy_cicle[Ninter,0]-energy_cicle[0,0])))
    if(plots==0):
        f4.show(); plt.show()

    # plots energies, density, and % in time
    if(plots==0):
        plot_convergence(tevol,energy_cicle,Ninter); plt.show()
        # other plots (for the final state)
        # plot_phase(z,psi,Zmax,t); plt.show()
        # plot_real_imag(z,psi,Zmax,t); plt.show()
        plot_wave_function(tevol, wave_function); plt.show()

    # closes files
    fn.close(); fe.close(); fpsi_all.close(); #fc.close()

    return c
