
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
# Program: First, evolution in imaginary time. We find an excited state of the condensate with soliton in the position that
# user decides. Second, evolution in real time. We evolve in real time the soliton, the user decides the number of oscilations.
# The program uses a local folder called "gpe_fft_utilities". There is different type of plots and an important function called
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
import numpy.linalg as lin
from pylab import*
import os, glob
try:
    from tkinter import *  
    from tkinter import ttk
    v3=True
except ImportError:
    from Tkinter import *
    import ttk
    v3=False
close('all')
pi=np.pi


# Data block
# ________________________________________________________________________________________________

# In[2]:

Zmax = 15.0              # Grid half length
Npoint =512              # Number of grid points
Nparticle = 500          # Number of particles
a_s = 0.5                # scattering length
whoz = 1.0               # harmonic oscilator angular frequency
Omega = pi/(2*Zmax)      # reference frame velocity
Ntime_fin = 1000         # total number of time steps
Ntime_out = 100          # number of time steps for intermediate outputs
Dtr = 1.0e-3             # real time step
Dti = 1.0e-3             # imaginary time step

# We choose the initial position of soliton:



class Demodark:
    def __init__(self,master,v3):
        self.name1=""
        self.name2=""
        self.master=master
        self.v3=v3
        frame=Frame(self.master)
        frame.pack()
        self.t_name1=DoubleVar()
        self.t_name2=DoubleVar()
        nb = ttk.Notebook(rootdark)
        page1 = ttk.Frame(nb)
        page2 = ttk.Frame(nb)        
        nb.add(page1, text='Programa')
        nb.add(page2, text='Notas')
        nb.pack(expand=1, fill="both")
        if self.v3==True:
            label1 = Label(page1, text="Posicion del soliton", background="black", foreground="white",font = "Verdana 16 bold")
            label1.pack(fill=X)
            bar1 = Scale(page1, from_=-6, to=6, variable=self.t_name1, length=600,tickinterval=2, resolution=0.1, orient=HORIZONTAL) 
            bar1.set(0)
            bar1.configure(bg='black',fg='white')
            bar1.pack(pady=10)
            label2 = Label(page1, text="Numero de oscilaciones", background="black", foreground="white",font = "Verdana 16 bold")
            label2.pack(fill=X)
            bar2 = Scale(page1, from_=1, to=10, variable=self.t_name2, length=600,tickinterval=3, resolution=1, orient=HORIZONTAL) 
            bar2.set(1)
            bar2.configure(bg='black',fg='white')
            bar2.pack()
        else:
            label1 = Label(page1, text="Posicion del soliton", background="black", foreground="white",font = "Verdana 16 bold")
            label1.pack(fill=X)
            bar1 = ttk.Scale(page1, from_=-6, to=6, variable=self.t_name1, length=600, orient=HORIZONTAL) 
            bar1.set(0)
            bar1.pack(pady=10)
            label2 = Label(page1, text="Numero de oscilaciones", background="black", foreground="white",font = "Verdana 16 bold")
            label2.pack(fill=X)
            bar2 = ttk.Scale(page1, from_=1, to=10, variable=self.t_name2, length=600, orient=HORIZONTAL) 
            bar2.set(1)
            bar2.pack()
        self.button=Button(page1, text='OK', command=self.show_values).pack()
        
        text2 = Text(page2, height=15, width=70)
  #      scroll = Scrollbar(root, command=text2.yview)
#        text2.configure(yscrollcommand=scroll.set)
        text2.tag_configure('bold_italics', font=('Arial', 12, 'bold', 'italic'))
        text2.tag_configure('big', font=('Verdana', 20, 'bold'))
        text2.tag_configure('color', foreground='#476042', 
						font=('Tempus Sans ITC', 12, 'bold'))
  #      text2.tag_bind('follow', '<1>', lambda e, t=text2: t.insert(END, "Not now, maybe later!"))
        text2.insert(END,'\nWilliam Shakespeare\n', 'big')
        quote = """
        To be, or not to be that is the question:
        Whether 'tis Nobler in the mind to suffer
        The Slings and Arrows of outrageous Fortune,
        Or to take Arms against a Sea of troubles,
        """
        text2.insert(END, quote, 'color')        
        text2.pack(side=LEFT)
        
   #     s = Separator(frame, orient=HORIZONTAL)
   #     s.set(0)


    def show_values(self):
        self.name1= (self.t_name1.get())
        self.name2= (self.t_name2.get())
        rootdark.destroy()
        
rootdark = Tk()
rootdark.wm_title("Dark solitons")
Dedark=Demodark(rootdark,v3)
rootdark.mainloop()


x0=Dedark.name1
osci=float(int(Dedark.name2))





#while True:
#    try:
#        x0=float(raw_input("introduce posicion inicial del soliton de 0 a 6"))
#        while (np.abs(x0)>(6)): # tolerance for the initial position of the soliton
#            print "ERROR: el soliton debe estar dentro del rango de 0 a 6"
#            x0=float(raw_input("introduce posicion inicial del soliton"))
#        break
#    except ValueError:
#        print("Escribe un numero")

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

def node(x,n,x0):# initial wave function
    fx = np.tanh(x-x0); # define the initial wf in R3
    return fft(fx)/n;   # FFT to K3

def normaliza(c): # normalization to 1
    norm = lin.norm(c)
    if ((norm-1.0)>1.0e-4): # check norm
        print("normalization from: ",norm)
    return c/norm


# Choose initial wave function and evolve in imaginary time:
# __________________________________________________________________________________________

# In[8]:

# initial wf: function defined at 'node', centered at x=x0
c0=normaliza(node(zp,Npoint,x0)); # wf at t=0
# evolve in time: parameters
t0=0.0
tevol=np.empty([Ninter+1])          # time vector
energy_cicle=np.empty([Ninter+1,5]) # put the energies in a matrix
energy_cicle[0,:] = Energy(c0)      # Energies at t=0
print("Energies in evolution imaginary time:          Emed    mu    Ekin    Epot    Eint")
print("         initial = %g %g %g %g %g"%(Energy(c0)))

#initiation of time evolution:
c=c0
tevol[0]=t0
j=0
t=0

#imaginary time evolution cicle:
for i in range(1, Ntime_fin+1):
    t += Dt.real
    psi=ifft(T_K(t0,-1j*Dt.real,c))*Npoint
    c=T_K(t0,-1j*Dt.real,fft(T_R_psi(t0,-1j*Dt.real,psi))/Npoint)
    c = normaliza(c); # check norm in the wf
    if(not(i%Ntime_out)):
        j+=1
        tevol[j] = t
        energy_cicle[j,:] = Energy(c)
print("         final = %g %g %g %g %g"%(Energy(c)))


# Plots of the final state:
#____________________________________________________________________________________________

# In[9]:

# Plot convergence during the evolution in the average energy per particle:

plot_convergence(tevol,energy_cicle[:,0],Ninter) # convergence in the average energy per particle

cc = ifft(c)*Npoint*NormWF**0.5      # FFT from K3 to R3 and include the wf norm
psi = changeFFTposition(cc,Npoint,0) # psi is the final wave function

psi*=np.exp(1j*pi/3) # This is useful to plot the wave function phase.
#plot different propieties of psi:

plot_density(z,psi,Zmax,t)
plot_phase(z,psi,Zmax,t)
plot_real_imag(z,psi,Zmax,t)


# New data block for real time evolution
#___________________________________________________________________________________________

# In[10]:

psi_sol=psi                  # we chance name variable
# User decide the oscilation number
#while True:
#    try:
#        osci =float(raw_input('introduce el numero de oscilaciones que quieres que de el soliton entre 1 y 10'))
#        while ((osci<1) or (osci>10)):
#            print "ERROR: las oscilaciones deben de estar entre un rango de 1 y 10"
#            osci=float(raw_input("introduce el numero de oscilaciones que quieres que de el soliton"))
#        break
#    except ValueError:
#        print("Escribe un numero")


Ntime_fin=int(osci*8800)     # total number of time steps
Ntime_out = 100              # number of time steps for intermediate outputs
Dtr=1.0e-3                   # real time step
Dti=1.0e-3                   # imaginary time step
Dt = Dtr-1j*Dti              # complex time
Ninter = Ntime_fin/Ntime_out # Number of outputs with the intermediate states

# Choose initial wave function and evolve in imaginary time:
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
tevol[0]=t0
j=0
t=0
tevol=np.empty([Ninter+1])     # time vector
pos_minus=np.empty([Ninter+1]) # put the minus position in a vector
val_minus=np.empty([Ninter+1]) # put the minus value in a vector
energi=np.empty([5])           # put the energies in a vector

# where the files of evolution will be saved
dir = "ds_evolution" # name of the directory
if (not os.path.exists(dir)): # creates the directory
    os.makedirs(dir)
else: # removes its contents if it already exists
    print ("Directory %r already exists. Do you want to continue?" % (dir))
    while True:
        try:
            if v3==True:
                ans = str(input("\tY/N: ")) # pause, answer to continue
            else:
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

file3=open('min.txt','w')
file3.write('Tabla donde se representan la posicion del minimo del soliton y su valor asociado en funcion del tiempo.\n')
file3.write('Tiempo\tPosicion del minimo\tValor del minimo\n')

#file2=open('WfDs.txt','w')
#file2.write('Datos de interes: N.particulas=%g\tPar.Interaccion=%g\tLong.caja=%g\tN.puntos=%g\tFreq.Oscilador=%g\tPot. quim.=%s\n ' %(Nparticle,gint,2*Zmax,Npoint,whoz,energi[1]))
#file2.write('x\tDensidad\tFase\tRe\tIm\tV(x)\n')

file4=open('phase.txt','w')

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
        file2=open('./%s/WfDs-%08d.txt'%(dir, j),'w')
        file2.write('Tiempo=%s\n' %(t))
        file2.write('Datos de interes: N.particulas=%g\tPar.Interaccion=%g\tLong.caja=%g\tN.puntos=%g\tFreq.Oscilador=%g\tPot. quim.=%s\n ' %(Nparticle,gint,2*Zmax,Npoint,whoz,energi[1]))
        file2.write('x\tDensidad\tFase\tRe\tIm\tV(x)\n')
        for i in range (0,int(2*Zmax/Dz)):
            file2.write("%s\t%s\t%s\t%s\t%s\t%s \n" %(z[i],(abs(psi)**2)[i],(np.angle(psi))[i],psi.real[i],psi.imag[i],changeFFTposition(Vpot_R,Npoint,0)[i]))
#        file2.write('\n\n')


# Minus of density (soliton)
        point= x0/Dz
        rang=np.empty([int((Npoint/2)+point+16)-int((Npoint/2)-point-15)])
        rang_2=np.empty([int((Npoint/2)+point+16)-int((Npoint/2)-point-15)])

        for i in range(len(z)):
            if (i>(int((Npoint/2)-point-16)) and i<(int((Npoint/2)+point+16))):
                rang[i-int((Npoint/2)-point-15)]=str((abs(psi)**2)[i])

        for i in range(len(rang)):
            if (min(rang)==rang[i]):
# Saves in a vector minus position/value and writes in a file. Also, writes the difference phase that creates soliton.
                j+=0
                pos_minus[j]=z[i+int((Npoint/2)-point-15)]
                dif_phase=np.angle(psi[i+int((Npoint/2)-point)])-np.angle(psi[i+int((Npoint/2)-point-30)])
                val_minus[j]=min(rang)
                file3.write('%s\t%s\t%s\n' %(t,pos_minus[j],val_minus[j]))
                file4.write('%s\t%s\n' %(t,dif_phase))

file.close()
file2.close()
file3.close()
file4.close()

# Prints final energy (soliton)
print("         final = %g %g %g %g %g"%(Energy(c)))

# Plots minus position and value
f5=plt.figure()
plt.title('Minus position/value',fontsize=15)
plt.xlabel('$t$',fontsize=15)
plt.xticks(np.arange(tevol[0],tevol[Ninter],int(tevol[Ninter]/5)))
plt.locator_params('y',nbins=3)
plt.plot(tevol,pos_minus, 'r--',label='$position-minus$') # plot minus position
plt.plot(tevol,val_minus*100.0, 'b--',label='$value-minus$') # plot minus value
plt.legend(fontsize=15)
f5.show()
plt.show()
