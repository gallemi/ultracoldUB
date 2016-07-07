
# coding: utf-8
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as lin

# Utilities
# -------------------------------------------------------------------------------------

# In[1]:

def normaliza(c): # normalization to 1
    norm = lin.norm(c)
#    if ((norm-1.0)>1.0e-4): # check norm
#        print("normalization from: ",norm)
    return c/norm


def normaliza_print(c): # normalization to 1
    norm = lin.norm(c)
    if ((norm-1.0)>1.0e-4): # check norm
        print("normalization from: ",norm**2)
    return c/norm  
#
# In[ ]:

def changeFFTposition(f,N,j):  # change the order in vectors from FFT
#    
# f(0...N-1) is the vector to order
# N is the vector dimension
# j is a switch indicating the change direction
# physical order is f=[(-(Zmax-Dz):Dz:-Dz) (0:Dz:Zmax) ]
# FFT order is f=[(0:Dk:kmax) (-(kmax-kz):kz:-kz)]
#
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


# Plot functions:
# -------------------------------------------------------------------------------------

# In[ ]:

def plot_convergence(x,y,n):
    f1=plt.figure()
    plt.title('Convergence',fontsize=15)
    plt.xlabel('time ($t \, \\omega_{ho}$)',fontsize=15)
    plt.ylabel('Energy per particle ($E/\\hbar \,\\omega_{ho}$)',fontsize=15)
    #plt.axis([-Zmax,Zmax,0, 8])
    plt.xticks(np.arange(0, x[n]+1,x[n]/5))
    plt.locator_params('y',nbins=3)
    plt.plot(x, y, 'r-')
    #plt.plot(z, psi, 'r.')
    f1.show()

# In[ ]:

def plot_density(z,psi,Lz,t):
    f2=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    #plt.ylabel('$\\psi\,(x)$',fontsize=15)
    #plt.axis([-Zmax,Zmax,0, 8])
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, abs(psi)**2, 'b--',label='$|\psi|^2$') # plot density
    plt.legend(fontsize=15)
    f2.show()

# In[ ]:

def plot_phase(z,psi,Lz,t):
    f3=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    #plt.ylabel('$\\psi\,(x)$',fontsize=15)
    #plt.axis([-Zmax,Zmax,0, 8])
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, np.angle(psi), 'b.',label='$Arg(\psi)$')
    plt.legend(fontsize=15)
    f3.show()
    
# In[ ]:

def plot_real_imag(z,psi,Lz,t):
    f3=plt.figure()
    plt.title('State at $t \,\\omega_{ho}=%g$'%(t),fontsize=15)
    plt.xlabel('$x/a_{ho}$',fontsize=15)
    #plt.ylabel('$\\psi\,(x)$',fontsize=15)
    #plt.axis([-Zmax,Zmax,0, 8])
    plt.xticks(np.arange(-Lz, Lz+1,Lz/2))
    plt.locator_params('y',nbins=3)
    plt.plot(z, psi.real, 'r.',label='real$(\psi)$')
    plt.plot(z, psi.imag, 'b--',label='imag$(\psi)$')
    plt.legend(fontsize=15)
    f3.show()    