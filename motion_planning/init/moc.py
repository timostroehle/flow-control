"""......................................................................................
bubble moc
......................................................................................"""
import math
import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

class Init:
    def __init__(self):
        self.K = 200
        self.method = "mp" # for midpoint-rule or "ie" for implicit euler
        self.problem = "inv"  #choose: "dir" or  "inv"
        """ loading data """
        self.path_to_data_init = "./init/data/data_init/"
        self.path_to_data_solution = "./init/data/data_solution/"
        self.a0 = np.loadtxt(self.path_to_data_init+"a0.dat")
        self.q0 = np.loadtxt(self.path_to_data_init+"q0.dat")
        self.ns = int(np.loadtxt(self.path_to_data_init+"param.dat")[0])
        self.nt = int(np.loadtxt(self.path_to_data_init+"param.dat")[1])
        self.c0 = np.loadtxt(self.path_to_data_init+"param.dat")[2]
        self.rho = np.loadtxt(self.path_to_data_init+"param.dat")[3]
        self.elspos = np.loadtxt(self.path_to_data_init+"elspos.dat", dtype="int")
        self.elsneg = np.loadtxt(self.path_to_data_init+"elsneg.dat", dtype="int")
        self.x_init = np.loadtxt(self.path_to_data_init+"xinit.dat")
        self.bc = np.loadtxt(self.path_to_data_init+"bc.dat",dtype="int").tolist()
        self.nkn = len(self.x_init)
        
class Bubbles:
    #init = Init().x_init
    def __init__(self,Y):
        self.Y = Y
        
    def p1(self):
        fhg = int(4*init.nkn)
        P1 = np.eye(fhg)[:,np.s_[[i for i in range(0,fhg) if i not in init.bc]]] 
        return P1

    def p2(self):
        fhg = int(4*init.nkn)
        P2 = np.eye(fhg)[:,np.s_[init.bc]]
        return P2

    def x_geg(self):
        x = self.Y.reshape(-1,1)
        xgeg = np.dot(self.p2().T,x)
        return xgeg

    def x_ges(self):
        x = self.Y.reshape(-1,1)
        xges = np.dot(self.p1().T,x)
        return xges

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

    def d_beta(self,s):
        """ gaussian distribution, maximum at s=0 """
        sig = 0.1
        s0 = 0.0
        db1 = -0.05*(1/10)*(2/sig**3)*(s-s0)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) 
        """ gaussian distribution of stiffness """
        sig = 0.25
        s0 = 0.5
        db = 0.5*(1/10)*(2/sig**3)*(s-s0)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) + db1
        """ constant distribution of stiffness """
        #db = 0
        return (1/(0+1))*db

    def dd_beta(self,s):
        """ gaussian distribution, maximum at s=0 """
        sig = 0.1
        s0 = 0.0
        a1 = 0.05*(1/10)
        ddb1 = -(2*a1/sig**3)*np.exp(-((s-s0)**2)/sig**2) + (4*a1*((s-s0)**2)/sig**5)*np.exp(-((s-s0)**2)/sig**2) 
        # """ gaussian distribution of stiffness """
        sig = 0.25
        s0 = 0.5 
        a = -0.5*(1/10)
        ddb = -(2*a/sig**3)*np.exp(-((s-s0)**2)/sig**2) + (4*a*((s-s0)**2)/sig**5)*np.exp(-((s-s0)**2)/sig**2) + ddb1
        """ constant distribution of stiffness """
        #ddb = 0
        return (1/(0+1))*ddb
    
    def f(self,xges):
        xgeg = self.x_geg()
        x = (np.dot(self.p1(),xges) + np.dot(self.p2(),xgeg)).reshape(-1,4)
        K = init.K
        f = []
        for i in range(0,len(init.elspos)):
            n = init.elspos[i,0]
            m = init.elspos[i,1]
            sm = x[m,0]; tm = x[m,1]; qm = x[m,2]; am = x[m,3]
            sn = x[n,0]; tn = x[n,1]; qn = x[n,2]; an = x[n,3]

            """ euler impl. """
            if init.method == "ie":
                b = self.beta(sm)
                db = self.d_beta(sm)
                g = am**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*am**(0.25))
                #dir
                f += [(sm-sn)/(tm-tn) - (qm/am +c)]
                #kom
                f += [(qm-qn)/(tm-tn) - (qm/am - c)*(am-an)/(tm-tn) + K*(qm/am) + a_m*(1/init.rho)*db*g]
            """ mid-point   """
            if init.method == "mp":
                s_mid = 0.5*(sm+sn)
                q_mid = 0.5*(qm+qn)
                a_mid = 0.5*(am+an)
                b = self.beta(s_mid)
                db = self.d_beta(s_mid)
                g = a_mid**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*a_mid**(0.25))
                #dir
                f += [(sm-sn)/(tm-tn)-(q_mid/a_mid+c)]
                #kom
                f += [(qm-qn)/(tm-tn)-(q_mid/a_mid-c)*(am-an)/(tm-tn)+K*(q_mid/a_mid) + a_mid*(1/init.rho)*db*g]
            else:
                print("numerical method is not implemented")
        for i in range(0,len(init.elsneg)):
            n = init.elsneg[i,0]
            m = init.elsneg[i,1]
            sm = x[m,0]; tm = x[m,1]; qm = x[m,2]; am = x[m,3]
            sn = x[n,0]; tn = x[n,1]; qn = x[n,2]; an = np.sqrt(x[n,3]**2)

            """ euler impl. """
            if init.method == "ie":
                b = self.beta(sm)
                db = self.d_beta(sm)
                g = am**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*am**(0.25))
                #dir
                f += [(sm-sn)/(tm-tn) - (qm/am - c)]
                #kom
                f += [(qm-qn)/(tm-tn) - (qm/am + c)*(am-an)/(tm-tn) + K*(qm/am) + a_m*(1/init.rho)*db*g]
            """ mid - point """
            if init.method == "mp":
                s_mid = 0.5*(sm+sn)
                q_mid = 0.5*(qm+qn)
                a_mid = 0.5*(am+an)
                b = self.beta(s_mid)
                db = self.d_beta(s_mid)
                g = a_mid**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*a_mid**(0.25))
                #dir
                f += [(sm-sn)/(tm-tn)-(q_mid/a_mid-c)]
                #kom
                f += [(qm-qn)/(tm-tn)-(q_mid/a_mid+c)*(am-an)/(tm-tn)+K*(q_mid/a_mid) + a_mid*(1/init.rho)*db*g]
            else:
                print("numerical method is not implemented")
        f = np.array(f).reshape(-1,1)
        return f

    def ana_tan(self,xges):
        xgeg = self.x_geg()
        x = (np.dot(self.p1(),xges) + np.dot(self.p2(),xgeg)).reshape(-1,4)
        K = init.K
        dr = np.zeros([1,4*init.nkn])
        for i in range(0,len(init.elspos)):
            n = init.elspos[i,0]
            m = init.elspos[i,1]
            sm = x[m,0]; tm = x[m,1]; qm = x[m,2]; am = x[m,3]
            sn = x[n,0]; tn = x[n,1]; qn = x[n,2]; an = x[n,3]
            """ euler impl. """
            if init.method == "ie":
                b = self.beta(sm)
                db = self.d_beta(sm)
                g = am**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*am**(0.25))
                #dir
                dr1 = np.zeros([1,4*init.nkn])
                dr = np.append(dr,dr1,axis=0)
                #kom
                dr1 = np.zeros([1,4*init.nkn])
                dr = np.append(dr,dr1,axis=0)
                try:
                    dr = jacobian_ie
                except ValueError:
                    print("Analytical Jacobian not implemented for Euler implicit...")
            """ mid-point """
            if init.method == "mp":
                s_mid = 0.5*(sm+sn)
                q_mid = 0.5*(qm+qn)
                a_mid = 0.5*(am+an)
                b = self.beta(s_mid)
                db = self.d_beta(s_mid)
                ddb = self.dd_beta(s_mid)
                g = a_mid**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*a_mid**(0.25))
                #dir
                dr1 = np.zeros([1,4*init.nkn])
                dr1[0,4*n+0] = -1/(tm-tn) - (init.rho**(-0.5))*0.25*(2**(-0.5))*(b**(-0.5))*db*(a_mid**(0.25))
                dr1[0,4*m+0] =  1/(tm-tn) - (init.rho**(-0.5))*0.25*(2**(-0.5))*(b**(-0.5))*db*(a_mid**(0.25))
                dr1[0,4*n+1] =  (sm-sn)/(tm-tn)**2
                dr1[0,4*m+1] = -(sm-sn)/(tm-tn)**2
                dr1[0,4*n+2] = -0.5*(1/a_mid)
                dr1[0,4*m+2] = -0.5*(1/a_mid)
                dr1[0,4*n+3] =  0.5*(q_mid/a_mid**2) - (init.rho**(-0.5))*(1/8)*(2**(-0.5))*(b**0.5)*a_mid**(-0.75)
                dr1[0,4*m+3] =  0.5*(q_mid/a_mid**2) - (init.rho**(-0.5))*(1/8)*(2**(-0.5))*(b**0.5)*a_mid**(-0.75)
                dr = np.append(dr,dr1,axis=0)
                #kom
                dr1 = np.zeros([1,4*init.nkn])
                dr1[0,4*n+0] = ((init.rho**(-0.5))*0.25*(2**(-0.5))*b**(-0.5)*db*a_mid**(0.25))*(am-an)/(tm-tn) + a_mid*(1/init.rho)*0.5*g*ddb 
                dr1[0,4*m+0] = ((init.rho**(-0.5))*0.25*(2**(-0.5))*b**(-0.5)*db*a_mid**(0.25))*(am-an)/(tm-tn) + a_mid*(1/init.rho)*0.5*g*ddb 
                dr1[0,4*n+1] =   (qm-qn)/(tm-tn)**2 - (q_mid/a_mid - c)*(am-an)/(tm-tn)**2 
                dr1[0,4*m+1] =  -(qm-qn)/(tm-tn)**2 + (q_mid/a_mid - c)*(am-an)/(tm-tn)**2
                dr1[0,4*n+2] = -1/(tm-tn) - 0.5*(1/a_mid)*(am-an)/(tm-tn) + 0.5*K*(1/a_mid)
                dr1[0,4*m+2] =  1/(tm-tn) - 0.5*(1/a_mid)*(am-an)/(tm-tn) + 0.5*K*(1/a_mid)
                dr1[0,4*n+3] =( (1/(tm-tn))*( (q_mid/a_mid) - c ) 
                               +(0.5*(q_mid/a_mid**2) + (init.rho**(-0.5))*0.125*(2**(-0.5))*(b**0.5)*(a_mid**(-0.75)))*(am-an)/(tm-tn)
                               - 0.5*K*(q_mid/a_mid**2)
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db)
                dr1[0,4*m+3] =(-(1/(tm-tn))*( (q_mid/a_mid) - c )
                               +(0.5*(q_mid/a_mid**2) + (init.rho**(-0.5))*0.125*(2**(-0.5))*(b**0.5)*(a_mid**(-0.75)))*(am-an)/(tm-tn)
                               - 0.5*K*(q_mid/a_mid**2)
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db)
                dr = np.append(dr,dr1,axis=0)
            else:
                print("numerical method is not implemented")
        for i in range(0,len(init.elsneg)):
            n = init.elsneg[i,0]
            m = init.elsneg[i,1]
            sm = x[m,0]; tm = x[m,1]; qm = x[m,2]; am = x[m,3]
            sn = x[n,0]; tn = x[n,1]; qn = x[n,2]; an = x[n,3]
            """ euler impl. """
            if init.method == "ie":
                b = self.beta(sm)
                db = self.d_beta(sm)
                g = am**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*am**(0.25))
                #dir
                dr1 = np.zeros([1,4*init.nkn])
                dr = np.append(dr,dr1,axis=0)
                #kom
                dr1 = np.zeros([1,4*init.nkn])
                dr = np.append(dr,dr1,axis=0)
                try:
                    dr = jacobian_ie
                except ValueError:
                    print("Analytical Jacobian not implemented for Euler implicit...")
            """ mid - point """
            if init.method == "mp":
                s_mid = 0.5*(sm+sn)
                q_mid = 0.5*(qm+qn)
                a_mid = 0.5*(am+an)
                b = self.beta(s_mid)
                db = self.d_beta(s_mid)
                ddb = self.dd_beta(s_mid)
                g = a_mid**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*a_mid**(0.25))
                #dir
                dr1 = np.zeros([1,4*init.nkn])
                dr1[0,4*n+0] = -1/(tm-tn) + (init.rho**(-0.5))*0.25*(2**(-0.5))*(b**(-0.5))*db*(a_mid**(0.25))
                dr1[0,4*m+0] =  1/(tm-tn) + (init.rho**(-0.5))*0.25*(2**(-0.5))*(b**(-0.5))*db*(a_mid**(0.25))
                dr1[0,4*n+1] =  (sm-sn)/(tm-tn)**2
                dr1[0,4*m+1] = -(sm-sn)/(tm-tn)**2
                dr1[0,4*n+2] = -0.5*(1/a_mid)
                dr1[0,4*m+2] = -0.5*(1/a_mid)
                dr1[0,4*n+3] =  0.5*(q_mid/a_mid**2) + (init.rho**(-0.5))*(1/8)*(2**(-0.5))*(b**0.5)*a_mid**(-0.75)
                dr1[0,4*m+3] =  0.5*(q_mid/a_mid**2) + (init.rho**(-0.5))*(1/8)*(2**(-0.5))*(b**0.5)*a_mid**(-0.75)
                dr = np.append(dr,dr1,axis=0)
                #kom
                dr1 = np.zeros([1,4*init.nkn])
                dr1[0,4*n+0] = (-(init.rho**(-0.5))*0.25*(2**(-0.5))*b**(-0.5)*db*a_mid**(0.25))*(am-an)/(tm-tn) + a_mid*(1/init.rho)*0.5*g*ddb
                dr1[0,4*m+0] = (-(init.rho**(-0.5))*0.25*(2**(-0.5))*b**(-0.5)*db*a_mid**(0.25))*(am-an)/(tm-tn) + a_mid*(1/init.rho)*0.5*g*ddb
                dr1[0,4*n+1] =   (qm-qn)/(tm-tn)**2 - (q_mid/a_mid + c)*(am-an)/(tm-tn)**2 
                dr1[0,4*m+1] =  -(qm-qn)/(tm-tn)**2 + (q_mid/a_mid + c)*(am-an)/(tm-tn)**2
                dr1[0,4*n+2] = -1/(tm-tn) - 0.5*(1/a_mid)*(am-an)/(tm-tn) + 0.5*K*(1/a_mid)
                dr1[0,4*m+2] =  1/(tm-tn) - 0.5*(1/a_mid)*(am-an)/(tm-tn) + 0.5*K*(1/a_mid)
                dr1[0,4*n+3] =( (1/(tm-tn))*( (q_mid/a_mid) + c )
                               +(0.5*(q_mid/a_mid**2) - (init.rho**(-0.5))*0.125*(2**(-0.5))*(b**0.5)*(a_mid**(-0.75)))*(am-an)/(tm-tn)
                               - 0.5*K*(q_mid/a_mid**2)
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db)
                dr1[0,4*m+3] =(-(1/(tm-tn))*( (q_mid/a_mid) + c )
                               +(0.5*(q_mid/a_mid**2) - (init.rho**(-0.5))*0.125*(2**(-0.5))*(b**0.5)*(a_mid**(-0.75)))*(am-an)/(tm-tn)
                               - 0.5*K*(q_mid/a_mid**2)
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db)
                dr = np.append(dr,dr1,axis=0)
            else:
                print("numerical method is not implemented")
        """ delete dirichlet bc """
        dr = np.delete(dr,0,0)
        #print(np.s_[i for i in init.bc])
        dr = np.delete(dr,np.s_[init.bc],1)
        return dr

    def num_tan(self,res,x):
        h_num =  1e-8
        dr = np.zeros([len(x),len(x)])
        for jj in range(0,len(x)):
            xh = np.array(x)
            xh[jj,0] = x[jj,0] + h_num
            res_h = self.f(xh)
            dr[:,jj] = (res_h[:,0]-res[:,0])/h_num
        return dr

    def solve(self,x):
        res_norm = 1
        while res_norm > 1e-8:
            print(res_norm)
            """ ..............................
            history = open('history.txt', 'a')
            history.write(str(res_norm))
            history.write('\n')
            history.close()
            ..............................."""
            """ resiudual vector .........."""
            res = self.f(x)
            """ jacobian .................."""
            dres = self.ana_tan(x)
            """ residual norm ............."""
            res_norm = np.linalg.norm(res)
            """ solve ....................."""
            dx = np.linalg.solve(dres,res)
            """ update ...................."""
            x = x - dx
        return x
    
""" initialisation ..................................................................."""
init = Init()
X0 = init.x_init
x = X0
X = np.array(x[0:init.ns+1,:])

bubb = Bubbles(x)
x0 = bubb.x_ges()

s = np.arange(0,1,0.01)
for i in range(1,len(s)):
    plt.plot([s[i-1],s[i]],[bubb.beta(s[i-1]),bubb.beta(s[i])],"-",color="b")
plt.grid()
plt.show()


""" Solve """
y = bubb.solve(x0)
    
x = (np.dot(bubb.p1(),y) + np.dot(bubb.p2(),bubb.x_geg())).reshape(-1,4)
""" plot initial characteristic net  """
for j in range(0,len(init.elspos[:,0])):
    mpos = init.elspos[j,0]
    npos = init.elspos[j,1]
    plt.plot([x[mpos,0],x[npos,0]],[x[mpos,1],x[npos,1]],"-",
                             color="tab:blue",linewidth=0.5)
for j in range(0,len(init.elsneg[:,0])):
    mneg = init.elsneg[j,0]
    nneg = init.elsneg[j,1]
    plt.plot([x[mneg,0],x[nneg,0]],[x[mneg,1],x[nneg,1]],"-",
                             color="tab:blue",linewidth=0.5)
plt.show()

""" xinit for characteristic step scheme  """
x0 = x[-(2*init.ns+1):,:]

""" save data """
np.savetxt(init.path_to_data_solution+"X.dat",x)
np.savetxt(init.path_to_data_solution+"x0.dat",x0)
