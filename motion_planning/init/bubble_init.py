"""......................................................................................
bubble init
......................................................................................"""
import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

class Init:
    def __init__(self):
        self.ns = int(sys.argv[1]) 
        self.nt = int(sys.argv[2])
        self.L = 1
        self.T = 1
        self.a0 = 0.1
        self.q0 = 0.0
        self.c0 = 5e4
        self.rho = 1.0e-6
        self.hs = self.L/self.ns
        
    def elspos(self):
        elspos = []
        k = 0
        """ init """
        for i in range(0,self.ns):
            k += 1
            for j in range(0,2*i+2):
                k+=1
                elspos += [k-1,k]
        """ steps """
        for i in range(0,self.nt-self.ns):
            k += 1
            for j in range(0,2*self.ns):
                k += 1
                elspos += [k-1,k]
        elspos = np.array(elspos).reshape(-1,2)
        return elspos

    def elsneg(self):
        elsneg = []
        k = 0
        """ init  """
        for i in range(0,self.ns):
            k += 1
            for j in range(0,2*i+1):
                k += 1
                elsneg += [k-(2+i),k+i]
        """ steps """
        for i in range(0,self.nt-self.ns):
            k += 1
            for j in range(0,2*self.ns):
                k+=1
                elsneg += [k-(self.ns+1),k+(self.ns-1)]
        elsneg = np.array(elsneg).reshape(-1,2)
        return elsneg

    def bc(self):
        bc = [0,1,2,3]
        k = 0
        """ init """
        for i in range(0,self.ns):
            k += 1
            bc += [4*k,4*k+1,4*k+2,4*k+3]
            for j in range(0,2*i+2):
                k+=1
            bc += [4*k,4*k+3]
        """ steps """
        for i in range(0,self.nt-self.ns):
            k += 1
            bc += [4*k]
            for j in range(0,2*self.ns):
                k += 1
            bc += [4*k,4*k+3]
        return bc

    def bc_trj(self):
        bc = []
        k = 0
        """ init """
        for i in range(0,self.ns):
            k += 1
            for j in range(0,2*i+2):
                k+=1
        """ steps """
        for i in range(0,self.nt-self.ns):
            k += 1
            for j in range(0,2*self.ns):
                k += 1
            bc += [k]
        return bc

    def beta(self,s):
        """ gaussian distribution, maximum at s=0 """
        sig = 0.1
        s0 = 0.0
        b1 = 0.05*(1/10)*(1/sig)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) +0.2
        """ gaussian distribution of stiffness """
        sig = 0.25
        s0 = 0.5
        b = -0.5*(1/10)*(1/sig)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) + 0.3 + b1 
        """ constant distribution of stiffness """
        #b = 2
        return (1/(0+1))*b + 0
    
    def c(self,s):
        b = self.beta(s)
        c = 0.0/self.a0 + 0.5*np.sqrt(2)*np.sqrt(self.beta(s))*self.a0**0.25
        """ const. velocity  """
        c = self.c0
        return c
    
    def xinit(self):
        xinit = [self.L, 0.0, self.q0, self.a0]
        """ init """
        for i in range(0,self.ns):
            t = 0.0
            s = self.L-(i+1)*self.hs
            ht = self.hs/self.c(s)
            xinit += [s, t , self.q0, self.a0]
            for j in range(1,2*i+3):
                s = self.L-(i+1)*self.hs+j*0.5*self.hs
                ht = self.hs/self.c(s)
                t += 0.5*ht
                xinit += [s, t , self.q0, self.a0]
        """ steps """
        for i in range(0,self.nt-self.ns):
            for j in range(0,2*self.ns+1):
                xinit += [j*0.5*self.hs, j*0.5*self.ht+(i+1)*self.ht, 0.0, self.a0]
        xinit = np.array(xinit).reshape(-1,4)
        return xinit

    def param(self):
        ns = self.ns
        nt = self.nt
        c0 = self.c0
        rho = self.rho
        param = np.array([ns,nt,c0,rho])
        return param
    
init = Init()

x = init.xinit()
""" plot initial characteristic net  """
for j in range(0,len(init.elspos()[:,0])):
    mpos = init.elspos()[j,0]
    npos = init.elspos()[j,1]
    plt.plot([x[mpos,0],x[npos,0]],[x[mpos,1],x[npos,1]],"-",
                             color="tab:blue",linewidth=0.5)
for j in range(0,len(init.elsneg()[:,0])):
    mneg = init.elsneg()[j,0]
    nneg = init.elsneg()[j,1]
    plt.plot([x[mneg,0],x[nneg,0]],[x[mneg,1],x[nneg,1]],"-",
                             color="tab:blue",linewidth=0.5)
plt.show()

np.savetxt("./init/data/data_init/elspos.dat",init.elspos(),fmt='%i')
np.savetxt("./init/data/data_init/elsneg.dat",init.elsneg(),fmt='%i')
np.savetxt("./init/data/data_init/xinit.dat",init.xinit())
np.savetxt("./init/data/data_init/bc.dat",init.bc(),fmt='%i')
np.savetxt("./init/data/data_init/bc_trj.dat",init.bc_trj(),fmt='%i')
np.savetxt("./init/data/data_init/param.dat",init.param())
np.savetxt("./init/data/data_init/a0.dat",np.array([init.a0]))
np.savetxt("./init/data/data_init/q0.dat",np.array([init.q0]))
