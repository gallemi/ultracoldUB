# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
from gpe_bs_parameters import *


def plot_convergence(x,y,n):
    f2=plt.figure()
    plt.title('Convergence',fontsize=15)
    plt.xlabel('time ($t \, \\omega_{\\xi}$)',fontsize=15)
    plt.ylabel('Energy per particle ($E/\\hbar \,\\omega_{\\xi}$)',fontsize=15)
    plt.xticks(np.arange(0, x[n]+1,x[n]/5))
    plt.locator_params('y',nbins=3)
    plt.plot(x, y[:,0], 'r-',label="$E_{med}$") # plot only average energy
    #plt.plot(x, y[:,2], 'b-',label="$E_{kin}$")
    #plt.plot(x, y[:,3], 'g-',label="$E_{pot}$")
    #plt.plot(x, y[:,4], 'y-',label="$E_{int}$")
    #plt.plot(x, y[:,1], 'm',label="$\\mu$")
    #plt.legend(fontsize=15)
    f2.show()

def plot_wave_function(x,y):
    f3=plt.figure()
    plt.title('Wave Function Integral',fontsize=15)
    plt.xlabel('time ($t \, \\omega_{\\xi}$)',fontsize=15)
    # plt.ylabel(' ',fontsize=15)
    plt.plot(x, y[:,0], 'r-', label='left side')
    plt.plot(x, y[:,1], 'g-', label='inside')
    plt.plot(x, y[:,2], 'b-', label='right side')
    plt.legend(fontsize=15)
    f3.show()

def plot_density(z,psi,Lz,t):
    f4=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, abs(psi)**2, 'b-',label='$|\psi|^2$') # plot density
    plt.legend(fontsize=15)
    f4.show()

def plot_phase(z,psi,Lz,t):
    f5=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, np.angle(psi), 'b.',label='$Arg(\psi)$')
    plt.legend(fontsize=15)
    f5.show()

def plot_real_imag(z,psi,Lz,t):
    f6=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, psi.real, 'r.',label='real$(\psi)$')
    plt.plot(z, psi.imag, 'b--',label='imag$(\psi)$')
    plt.legend(fontsize=15)
    f6.show()
