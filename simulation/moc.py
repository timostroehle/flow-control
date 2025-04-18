"""......................................................................................
bubble moc
......................................................................................"""
import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# from sympy import *
# init_printing()
np.set_printoptions(precision=3, threshold=1e5, linewidth=300)

class Init:
    def __init__(self):
        self.ns = int(sys.argv[1])  # choose even number (-> centered bubble)
        self.L = 1.0
        self.T = 0.05
        self.c0 = 1000#165.5# 165.5#0.2#1#1.0#1.8 #1.4
        self.hs = self.L/self.ns
        self.ht = self.hs/self.c0
        self.nt = 1#1 #int(np.ceil(self.T/self.ht))
        self.nel = 2*self.ns*self.nt
        self.nkn = (self.ns*self.nt)+(self.ns+1)*(self.nt+1)
        self.K = 100#100#25#1.00
        self.rho = 1e-6
        self.a0 = 0.1#1.0 #1.0
        self.da = 0.001#0.1
        self.q0 = 0.0
        self.dq = 0.0# 0.01#0.1
        self.method = "mp" # "mp" # mp for midpoint-rule; alternatively "ie" for implicit euler
        self.integration = "local"  #choose: "local" or  "global"

    def elspos(self):
        els = np.zeros([self.nel,2],dtype=int)
        k = 0
        for j in range(0,self.nt):
            for i in range(0,self.ns):
                els[k,:] = [i+j*(2*(self.ns+1)-1),i+j*(2*(self.ns+1)-1)+self.ns+1]
                els[k+1,:] = [i+j*(2*(self.ns+1)-1)+self.ns+1,i+j*(2*(self.ns+1)-1)+2*self.ns+2]
                k += 2
        return els

    def elsneg(self):
        els = np.zeros([self.nel,2],dtype=int)
        k = 0
        for j in range(0,self.nt):
            for i in range(0,self.ns):
                els[k,:] = [i+j*(2*(self.ns+1)-1)+1,i+j*(2*(self.ns+1)-1)+self.ns+1]
                els[k+1,:] = [i+j*(2*(self.ns+1)-1)+self.ns+1,i+j*(2*(self.ns+1)-1)+2*self.ns+2-1]
                k += 2 
        return els

    def st_init(self):
        st = np.array([[0,0]])
        for j in range(0,self.nt):
            for i in range(0,self.ns+1):
                st = np.append(st,np.array([[0+i*self.hs , j*self.ht]]),axis=0)
            for i in range(0,self.ns):
                st = np.append(st,np.array([[self.hs*(i+0.5), (j+0.5)*self.ht]]),axis=0)
        for i in range(0,self.ns+1):
            st = np.append(st,np.array([[0+i*self.hs, self.nt*self.ht]]),axis=0)
        st = np.delete(st,0,0)
        return st
    
    def x_init(self):
        st = self.st_init()
        q = 0*np.ones([len(st),1])
        a = self.a0*np.ones([len(st),1])
        a[2*self.ns+1,0] = self.a_pre(st[2*self.ns+1,1])
        X0 = np.append(np.append(st,q,axis=1),a,axis=1)
        return X0

    def bc_init(self):
        X0 = self.x_init()
        bc = []
        for i in range(0,len(X0)):
            if X0[i,1]==0:
                bc += [4*i, 4*i+1, 4*i+2, 4*i+3]
            elif X0[i,0]==0 and X0[i,1]!=0:
                bc += [4*i, 4*i+3]  #"""direct imposition of dirichlet bc for a on s=0,t"""
                #bc += [4*i]
            elif X0[i,0]==1 and X0[i,1]!=0:
                bc += [4*i, 4*i+3]
        return bc

    def a_pre(self,t):
        t0 = 0.00; t1 = 0.020; dt = t1-t0; da = self.da #t1=0.14
        if t<t0:
            a = self.a0
        elif t>t1:
            a = self.a0+da
        else:
            """ C2-function (cos(cos())) """
            f = 0.5* (np.cos((t-t0)*np.pi/dt + np.pi) + 1)
            a = da*(0.5*np.cos(f*np.pi+np.pi) + 0.5) + 1
            """ C0 linear """
            #p = (1/(tf-t0))*t - t0/(tf-t0)
            """ C-1 constant """
            #a = init.a0 + da
        """ sigmoid """
        # a = da*(1+np.exp(-1*(t-10)))**(-1) + 1
        return a

class Bubbles:
    init = Init().x_init()
    def __init__(self,Y,i):
        self.Y = Y
        self.i = i
        
    def p1(self):
        fhg = int(4*init.nkn)
        P1 = np.eye(fhg)[:,np.s_[[i for i in range(0,fhg) if i not in init.bc_init()]]] 
        return P1

    def p2(self):
        fhg = int(4*init.nkn)
        P2 = np.eye(fhg)[:,np.s_[init.bc_init()]]
        return P2

    # def st(self):
    #     st = np.array([[0,0]])
    #     for j in range(0,init.nt):
    #         for i in range(0,init.ns+1):
    #             st = np.append(st,np.array([[0+i*init.hs , (j+self.i)*init.ht]]),axis=0)
    #         for i in range(0,init.ns):
    #             st = np.append(st,np.array([[init.hs*(i+0.5), (j+0.5+self.i)*init.ht]]),axis=0)
    #     for i in range(0,init.ns+1):
    #         st = np.append(st,np.array([[0+i*init.hs , (init.nt+self.i)*init.ht]]),axis=0)
    #     st = np.delete(st,0,0)
    #     return st

    
    # """direct imposition of dirichlet bc for a on s=0,t"""
    def X(self):
        X = self.Y
        x[2*init.ns+1,3] = self.a_pre(x[2*init.ns+1,1])
        return X
    # """................................................"""

    def x_geg(self):
        x = self.X().reshape(-1,1) ### direct imposition of dirichlet bc for a on s=0,t"""
        #x = self.Y.reshape(-1,1)
        xgeg = np.dot(self.p2().T,x)
        return xgeg

    def x_ges(self):
        #x = self.X().reshape(-1,1) """direct imposition of dirichlet bc for a on s=0,t"""
        x = self.Y.reshape(-1,1)
        xges = np.dot(self.p1().T,x)
        return xges

    # def c_guess(self,s,q,a):
    #     b = self.beta(s)
    #     c = q/a - (init.rho**(-0.5))*0.5*np.sqrt(2)*np.sqrt(self.beta(s))*a**0.25
    #     """ alternative: c = const. in s  """
    #     #c = self.c0
    #     return c
    
    def update(self,X):
        A = np.array(X)
        for i in range(0,init.ns+1):
            A[i,0] = X[i+2*init.ns+1,0]
            A[i,1] = X[i+2*init.ns+1,1]
            A[i,2] = X[i+2*init.ns+1,2]
            A[i,3] = X[i+2*init.ns+1,3]

        for i in range(0,init.ns):
            hti = X[i+2*init.ns+1,1] - X[i,1]
            A[i+init.ns+1,0]  = X[i+init.ns+1,0]
            A[i+init.ns+1,1]  = X[i+init.ns+1,1] + hti #+ init.ht
            A[i+init.ns+1,2]  = 0.5*(X[i+2*init.ns+1,2] + X[i+2*init.ns+2,2])
            A[i+init.ns+1,3]  = 0.5*(X[i+2*init.ns+1,3] + X[i+2*init.ns+2,3])

            # A[i+self.ns+1,0]  = X[i+self.ns+1,0]
            # A[i+self.ns+1,1]  = X[i+self.ns+1,1] + self.ht
            # A[i+self.ns+1,2]  = X[i+2*self.ns+1,2]
            # A[i+self.ns+1,3]  = X[i+2*self.ns+1,3]
            
        for i in range(0,init.ns+1):
            hti = X[i+2*init.ns+1,1] - X[i,1]
            A[i+2*init.ns+1,0] = X[i+2*init.ns+1,0]
            A[i+2*init.ns+1,1] = X[i+2*init.ns+1,1] + hti #+ init.ht #+hti #
            A[i+2*init.ns+1,2] = X[i+2*init.ns+1,2] 
            A[i+2*init.ns+1,3] = X[i+2*init.ns+1,3]
        A[2*init.ns+1,3] = self.a_pre(X[2*init.ns+1,1])
        return A

    def a_pre(self,t):
        t0 = 0.0; t1 = 0.020; dt = t1-t0; da = init.da
        if t<t0:
            a = init.a0
        elif t>t1:
            a = init.a0+da
        else:
            """ C2-function (cos(cos())) """
            f = 0.5* (np.cos((t-t0)*np.pi/dt + np.pi) + 1)
            a = da*(0.5*np.cos(f*np.pi+np.pi) + 0.5) + init.a0
            """ C0 linear """
            #p = (1/(tf-t0))*t - t0/(tf-t0)
            """ C-1 constant """
            #a = init.a0 + da
        """ sigmoid """
        # a = da*(1+np.exp(-1*(t-10)))**(-1) + 1
        return a

    def da_pre(self,t):
        t0 = 0.0; t1 = 0.020; dt = t1-t0; da = init.da
        if t<t0:
            dq = 0
        elif t>t1:
            dq = 0
        else:
            """ C2-function (cos(cos())) """
            f = 0.5* (np.cos((t-t0)*np.pi/dt + np.pi) + 1)
            df = -0.5*np.sin((t-t0)*np.pi/dt + np.pi) * np.pi/dt
            dq = -0.5*da*np.sin(f*np.pi + np.pi)*np.pi*df 
            """ C0 linear """
            #dq = a*(1/(t1-t0))
            """ C-1 constant step"""
            #dq = 0
        """ sigmoid """
        # dq = ...
        """ load from data """
        return dq

    # def q_pre(self,t):
    #     t0 = 2.0; t1 = 5.0; dt = t1-t0; dq = init.dq
    #     if t<t0:
    #         q = init.q0
    #     elif t>t1:
    #         q = init.q0+dq
    #     else:
    #         """ C2-function (cos(cos())) """
    #         f = 0.5* (np.cos((t-t0)*np.pi/dt + np.pi) + 1)
    #         q = dq*(0.5*np.cos(f*np.pi+np.pi) + 0.5) 
    #         """ C0 linear """
    #         #p = (1/(tf-t0))*t - t0/(tf-t0)
    #         """ C-1 constant """
    #         #a = init.a0 + da
    #     """ sigmoid """
    #     # a = da*(1+np.exp(-1*(t-10)))**(-1) + 1
    #     return q

    def beta(self,s):
        """ gaussian distribution, maximum at s=0 """
        sig = 0.1 #5#*np.sqrt(2*np.pi)
        s0 = 0.0
        a1 = 0.05*(1/10)
        b1 = a1*(1/sig)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) + 0.2 #+ 2.0
        """ gaussian distribution of stiffness """
        sig = 0.25 #5#*np.sqrt(2*np.pi)
        s0 = 0.5
        a = -0.50*(1/10)
        b = a*(1/sig)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) + 0.3  + b1  # 2.5
        """ constant distribution of stiffness """
        #b = 0.1#0.05
        return b

    def d_beta(self,s):
        """ gaussian distribution, maximum at s=0 """
        sig = 0.1 #5#*np.sqrt(2*np.pi)
        s0 = 0.0
        a1 = 0.05*(1/10)
        db1 = -a1*(2/sig**3)*(s-s0)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) 
        """ gaussian distribution of stiffness """
        sig = 0.25#5#*np.sqrt(2*np.pi)
        s0 = 0.5
        a = -0.50*(1/10)
        db = -a*(2/sig**3)*(s-s0)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) + db1
        """ constant distribution of stiffness """
        #db = 0
        return db

    def dd_beta(self,s):
        """ gaussian distribution, maximum at s=0 """
        sig = 0.1#1.0 #0.1 #5#*np.sqrt(2*np.pi)
        s0 = 0.0
        a1 = 0.05*(1/10)
        ddb1 = -(2*a1/sig**3)*np.exp(-((s-s0)**2)/sig**2) + (4*a1*((s-s0)**2)/sig**5)*np.exp(-((s-s0)**2)/sig**2)
        # """ gaussian distribution of stiffness """
        sig = 0.25#2.0 #0.25#5#*np.sqrt(2*np.pi)
        s0 = 0.5#5.0 #0.5
        a = -0.5*(1/10)#2
        ddb = -(2*a/sig**3)*np.exp(-((s-s0)**2)/sig**2) + (4*a*((s-s0)**2)/sig**5)*np.exp(-((s-s0)**2)/sig**2) + ddb1
        """ constant distribution of stiffness """
        #ddb = 0
        return ddb
        
    def f(self,xges):
        xgeg = self.x_geg()
        x = (np.dot(self.p1(),xges) + np.dot(self.p2(),xgeg)).reshape(-1,4)
        K = init.K
        f = []
        for i in range(0,len(init.elspos())):
            n = init.elspos()[i,0]
            m = init.elspos()[i,1]
            sm = x[m,0]; tm = x[m,1]; qm = x[m,2]; am = x[m,3]
            sn = x[n,0]; tn = x[n,1]; qn = x[n,2]; an = x[n,3]
            #euler impl.
            if init.method == "ie":
                b = self.beta(sm)
                db = self.d_beta(sm)
                g = am**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*am**(0.25))
                #dir
                f += [(sm-sn)/(tm-tn) - (qm/am +c)]
                #kom
                f += [(qm-qn)/(tm-tn) - (qm/am - c)*(am-an)/(tm-tn) + K*(qm/am) + am*db*g*(1/init.rho)]
            #mid-point
            elif init.method == "mp":
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
                f += [(qm-qn)/(tm-tn)-(q_mid/a_mid-c)*(am-an)/(tm-tn)+K*(q_mid/a_mid) + a_mid*db*g*(1/init.rho)]
            else:
                print("numerical method is not implemented")
        for i in range(0,len(init.elsneg())):
            n = init.elsneg()[i,0]
            m = init.elsneg()[i,1]
            sm = x[m,0]; tm = x[m,1]; qm = x[m,2]; am = x[m,3]
            sn = x[n,0]; tn = x[n,1]; qn = x[n,2]; an = x[n,3]
            #euler impl.
            if init.method == "ie":
                b = self.beta(sm)
                db = self.d_beta(sm)
                g = am**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*am**(0.25))
                #dir
                f += [(sm-sn)/(tm-tn) - (qm/am - c)]
                #kom
                f += [(qm-qn)/(tm-tn) - (qm/am + c)*(am-an)/(tm-tn) + K*(qm/am) + am*db*g*(1/init.rho)]
            #mid - point
            elif init.method == "mp":
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
                f += [(qm-qn)/(tm-tn)-(q_mid/a_mid+c)*(am-an)/(tm-tn)+K*(q_mid/a_mid) + a_mid*db*g*(1/init.rho)]
            else:
                print("numerical method is not implemented")
            """ ........... """
        """indirect imposition of dirichlet bc for a on s=0,t"""
        # for i in range(1,init.nt+1):
        #     print(self.a_pre(x[i*(2*init.ns+1),1]))
        #     f += [x[i*(2*init.ns+1),3] - self.a_pre(x[i*(2*init.ns+1),1])]
            
        f = np.array(f).reshape(-1,1)
        return f

    def ana_tan(self,xges):
        xgeg = self.x_geg()
        x = (np.dot(self.p1(),xges) + np.dot(self.p2(),xgeg)).reshape(-1,4)
        K = init.K
        dr = np.zeros([1,4*init.nkn])
        for i in range(0,len(init.elspos())):
            n = init.elspos()[i,0]
            m = init.elspos()[i,1]
            sm = x[m,0]; tm = x[m,1]; qm = x[m,2]; am = x[m,3]
            sn = x[n,0]; tn = x[n,1]; qn = x[n,2]; an = x[n,3]
            #euler impl.
            if init.method == "ie":
                b = self.beta(sm)
                db = self.d_beta(sm)
                ddb = self.dd_beta(sm)
                g = am**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*am**(0.25))
                #dir
                dr1 = np.zeros([1,4*init.nkn])
                dr1[0,4*n+0] = -1/(tm-tn) #- (init.rho**(-0.5))*0.5*(2**(-0.5))*(b**(-0.5))*db*(am**(0.25))
                dr1[0,4*m+0] =  1/(tm-tn) - (init.rho**(-0.5))*0.5*(2**(-0.5))*(b**(-0.5))*db*(am**(0.25))
                dr1[0,4*n+1] =  (sm-sn)/(tm-tn)**2
                dr1[0,4*m+1] = -(sm-sn)/(tm-tn)**2
                dr1[0,4*n+2] =  0#-(1/am)
                dr1[0,4*m+2] = -(1/am)
                dr1[0,4*n+3] =  0#(qm/am**2) - (init.rho**(-0.5))*(1/4)*(2**(-0.5))*(b**0.5)*am**(-0.75)
                dr1[0,4*m+3] =  (qm/am**2) - (init.rho**(-0.5))*(1/4)*(2**(-0.5))*(b**0.5)*am**(-0.75)
                dr = np.append(dr,dr1,axis=0)
                #kom
                dr1 = np.zeros([1,4*init.nkn])
                dr1[0,4*n+0] = 0#((init.rho**(-0.5))*0.5*(2**(-0.5))*b**(-0.5)*db*am**(0.25))*(am-an)/(tm-tn)+ am*(1/init.rho)*g*ddb 
                dr1[0,4*m+0] = ((init.rho**(-0.5))*0.5*(2**(-0.5))*b**(-0.5)*db*am**(0.25))*(am-an)/(tm-tn)+ am*(1/init.rho)*g*ddb 
                dr1[0,4*n+1] =   (qm-qn)/(tm-tn)**2 - (qm/am - c)*(am-an)/(tm-tn)**2 
                dr1[0,4*m+1] =  -(qm-qn)/(tm-tn)**2 + (qm/am - c)*(am-an)/(tm-tn)**2
                dr1[0,4*n+2] = -1/(tm-tn) #- (1/am)*(am-an)/(tm-tn) + K*(1/am)
                dr1[0,4*m+2] =  1/(tm-tn) - (1/am)*(am-an)/(tm-tn) + K*(1/am)
                dr1[0,4*n+3] =( (1/(tm-tn))*( (qm/am) - c ) 
                               #+((qm/am**2) + (init.rho**(-0.5))*0.25*(2**(-0.5))*(b**0.5)*(am**(-0.75)))*(am-an)/(tm-tn)
                               #- K*(qm/am**2)
                               #+ (1/init.rho)*0.5*am**(0.5)*db + (1/init.rho)*g*db
                               )
                dr1[0,4*m+3] =(-(1/(tm-tn))*( (qm/am) - c )
                               +((qm/am**2) + (init.rho**(-0.5))*0.25*(2**(-0.5))*(b**0.5)*(am**(-0.75)))*(am-an)/(tm-tn)
                               - K*(qm/am**2)
                               + (1/init.rho)*0.5*am**(0.5)*db + (1/init.rho)*g*db
                               )
                dr = np.append(dr,dr1,axis=0)
            #mid-point 
            elif init.method == "mp":
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
                dr1[0,4*n+0] = ((init.rho**(-0.5))*0.25*(2**(-0.5))*b**(-0.5)*db*a_mid**(0.25))*(am-an)/(tm-tn)+ a_mid*(1/init.rho)*0.5*g*ddb 
                dr1[0,4*m+0] = ((init.rho**(-0.5))*0.25*(2**(-0.5))*b**(-0.5)*db*a_mid**(0.25))*(am-an)/(tm-tn)+ a_mid*(1/init.rho)*0.5*g*ddb 
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
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db )
                dr = np.append(dr,dr1,axis=0)
            else:
                print("numerical method is not implemented")
        for i in range(0,len(init.elsneg())):
            n = init.elsneg()[i,0]
            m = init.elsneg()[i,1]
            sm = x[m,0]; tm = x[m,1]; qm = x[m,2]; am = x[m,3]
            sn = x[n,0]; tn = x[n,1]; qn = x[n,2]; an = x[n,3]
            #euler impl.
            if init.method == "ie":
                b = self.beta(sm)
                db = self.d_beta(sm)
                ddb = self.dd_beta(sm)
                g = am**0.5 - init.a0**0.5
                c = (init.rho**(-0.5))*(2**(-0.5))*((b**0.5)*am**(0.25))
                #dir
                dr1 = np.zeros([1,4*init.nkn])
                dr1[0,4*n+0] = -1/(tm-tn) #+ (init.rho**(-0.5))*0.5*(2**(-0.5))*(b**(-0.5))*db*(am**(0.25))
                dr1[0,4*m+0] =  1/(tm-tn) + (init.rho**(-0.5))*0.5*(2**(-0.5))*(b**(-0.5))*db*(am**(0.25))
                dr1[0,4*n+1] =  (sm-sn)/(tm-tn)**2
                dr1[0,4*m+1] = -(sm-sn)/(tm-tn)**2
                dr1[0,4*n+2] =  0#-(1/am)
                dr1[0,4*m+2] = -(1/am)
                dr1[0,4*n+3] =  0#(qm/am**2) + (init.rho**(-0.5))*(1/4)*(2**(-0.5))*(b**0.5)*am**(-0.75)
                dr1[0,4*m+3] =  (qm/am**2) + (init.rho**(-0.5))*(1/4)*(2**(-0.5))*(b**0.5)*am**(-0.75)
                dr = np.append(dr,dr1,axis=0)
                #kom
                dr1 = np.zeros([1,4*init.nkn])
                dr1[0,4*n+0] = 0#(-(init.rho**(-0.5))*0.5*(2**(-0.5))*b**(-0.5)*db*am**(0.25))*(am-an)/(tm-tn) + am*(1/init.rho)*g*ddb
                dr1[0,4*m+0] = (-(init.rho**(-0.5))*0.5*(2**(-0.5))*b**(-0.5)*db*am**(0.25))*(am-an)/(tm-tn) + am*(1/init.rho)*g*ddb
                dr1[0,4*n+1] =   (qm-qn)/(tm-tn)**2 - (qm/am + c)*(am-an)/(tm-tn)**2 
                dr1[0,4*m+1] =  -(qm-qn)/(tm-tn)**2 + (qm/am + c)*(am-an)/(tm-tn)**2
                dr1[0,4*n+2] = -1/(tm-tn) #- (1/am)*(am-an)/(tm-tn) + K*(1/am)
                dr1[0,4*m+2] =  1/(tm-tn) - (1/am)*(am-an)/(tm-tn) + K*(1/am)
                dr1[0,4*n+3] =( (1/(tm-tn))*( (qm/am) + c )
                               #+((qm/am**2) - (init.rho**(-0.5))*0.25*(2**(-0.5))*(b**0.5)*(am**(-0.75)))*(am-an)/(tm-tn)
                               #- K*(qm/am**2)
                               #+ (1/init.rho)*0.5*am**0.5*db + (1/init.rho)*g*db
                               )
                dr1[0,4*m+3] =(-(1/(tm-tn))*( (qm/am) + c )
                               +((qm/am**2) - (init.rho**(-0.5))*0.25*(2**(-0.5))*(b**0.5)*(am**(-0.75)))*(am-an)/(tm-tn)
                               - K*(qm/am**2)
                               + (1/init.rho)*0.5*am**0.5*db + (1/init.rho)*g*db
                               )
                dr = np.append(dr,dr1,axis=0)
            #mid - point
            elif init.method == "mp":
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
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db )
                dr1[0,4*m+3] =(-(1/(tm-tn))*( (q_mid/a_mid) + c )
                               +(0.5*(q_mid/a_mid**2) - (init.rho**(-0.5))*0.125*(2**(-0.5))*(b**0.5)*(a_mid**(-0.75)))*(am-an)/(tm-tn)
                               - 0.5*K*(q_mid/a_mid**2)
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db )
                dr = np.append(dr,dr1,axis=0)
            else:
                print("numerical method is not implemented")
        """ algebraic constraints """
        # dr1 = np.zeros([1,4*init.nkn])
        # dr1[0,4*(2*init.ns+1)+3] = 1.0
        # dr1[0,4*(2*init.ns+1)+1] = -self.da_pre(x[2*init.ns+1,1])
        #dr = np.append(dr,dr1,axis=0)
        """ delete dirichlet bc """
        dr = np.delete(dr,0,0)
        #print(np.s_[i for i in init.bc_init()])
        dr = np.delete(dr,np.s_[init.bc_init()],1)
        return dr

    # def sym_tan(self):
    #     dr = 0
    #     return dr

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
        while res_norm > 1e-5:
            print(res_norm)
            """ ..............................
            history = open('history.txt', 'a')
            history.write(str(res_norm))
            history.write('\n')
            history.close()
            ..............................."""
            res = self.f(x)
            #dres = self.num_tan(res,x)
            dres = self.ana_tan(x)
            #dres = self.sym_tan()
            res_norm = np.linalg.norm(res)
            ### Normal:
            dx = np.dot(np.linalg.inv(dres),res)
            ### Thikonov:
            # G = 1e-12*np.eye(len(dr1))
            # AA = np.dot( np.linalg.inv( np.dot(np.transpose(dr1),dr1) + G) , np.transpose(dr1) )
            # dx1 = np.dot(AA,r1)
            x = x - dx
        return x

    def res_q_steady_state(self,q):
        K = init.K
        a1 = 1.1
        a0 = init.a0
        a2 = a0
        p1 = a1**0.5 - a0**0.5
        p2 = a2**0.5 - a0**0.5
        da_ds_1 = -(K*q/a1)/(0.5*3*a1**(0.5) - q**2/a1**2)
        da_ds_2 = -(K*q/a2)/(0.5*3*a2**(0.5) - q**2/a2**2)
        F1 = 0.25/a1**2 #np.log(a1)*(da_ds_1)**(-1)
        F2 = 0.25/a2**2 #np.log(a2)*(da_ds_2)**(-1)

        #R = K*q*(F2-F1)
        R = 0.25*K*((q/a2**2)+(q/a1**2))
        #R = -(1/q)*(p1-p2)
        res = 0.5*((q/a2)**2-(q/a1)**2) - (p1-p2) + R  
        return res
    
    def f_steady_state(self,s,a):
        K = init.K
        q0 = 1.0
        q = scipy.optimize.newton(self.res_q_steady_state,q0)
        z = (K*q/a + a*self.d_beta(s)*(a**0.5 - init.a0**0.5)  )
        n = 0.5*self.beta(s)*a**(0.5) - q**2/a**2
        f = -z/n
        return f

    def steady_state(self,a):
        sol = scipy.integrate.solve_ivp(self.f_steady_state,(0,1),[1.1],method="BDF",t_eval=s)
        a = sol.y[0]
        return a

    # def integrate_c(self,q):
    #     q = Symbol('q')
    #     a = Symbol('a')
    #     f = q/a
    #     g = 0.5*a*a**(-0.5) - q**2/a**2 
    #     c = 0.5*f + (0.25*f**2 + g)**0.5
    #     C = integrate(c, q)
    #     return C

    def riemann_pos(self,x):
        h = np.zeros([len(x),3])
        for i in range(0,len(h)):
            q = x[i,2]
            a = x[i,3]    
            h[i,0] = x[i,0]
            h[i,1] = x[i,1]
            h[i,2] = a - a*q + 0.5*q**2/a - 2**(-0.5)*q*a**(0.25) + 1
        return h

""" initialisation ..................................................................."""
init = Init()
X0 = init.x_init()
x = X0
X = np.array(x[0:init.ns+1,:])
G = np.zeros([1,2])

""" plot initial characteristic net  """
for j in range(0,len(init.elspos()[:,0])):
    mpos = init.elspos()[j,0]
    npos = init.elspos()[j,1]
    mneg = init.elsneg()[j,0]
    nneg = init.elsneg()[j,1]
    plt.plot([x[mpos,0],x[npos,0]],[x[mpos,1],x[npos,1]],"-",
                             color="tab:blue",linewidth=0.5)
    plt.plot([x[mneg,0],x[nneg,0]],[x[mneg,1],x[nneg,1]],"-",
                             color="tab:blue",linewidth=0.5)
plt.show()
""" plot beta-law """
bubb = Bubbles(x,0)
s = np.arange(0,1,0.01)
for i in range(1,len(s)):
    plt.plot([s[i-1],s[i]],[bubb.beta(s[i-1]),bubb.beta(s[i])],"-",color="b")
    #plt.plot([s[i-1],s[i]],[bubb.d_beta(s[i-1]),bubb.d_beta(s[i])],"-",color="r")
    #plt.plot([s[i-1],s[i]],[bubb.dd_beta(s[i-1]),bubb.dd_beta(s[i])],"-",color="g")
plt.grid()
plt.show()

s = np.arange(0,init.T,0.001)
for i in range(1,len(s)):
    plt.plot([s[i-1],s[i]],[bubb.a_pre(s[i-1]),bubb.a_pre(s[i])],"-",color="b")
    #plt.plot([s[i-1],s[i]],[bubb.d_beta(s[i-1]),bubb.d_beta(s[i])],"-",color="r")
    #plt.plot([s[i-1],s[i]],[bubb.dd_beta(s[i-1]),bubb.dd_beta(s[i])],"-",color="g")
plt.grid()
plt.show()

if init.integration == "global":
    tf = 1e-3 # assuring tf<T
else:
    tf = x[init.ns+(2*init.ns+1)*init.nt,1]

i = 0
while tf<init.T:
    # bubb = Bubbles(x,i)
    # s = [i*0.01 for i in range(0,100)]
    # a = bubb.steady_state(s)

    print("time-slab no." + str(i) + "; tf = " + str(tf))
    bubb = Bubbles(x,i)
    x0 = bubb.x_ges()
    """ Solve """
    #y = scipy.optimize.newton_krylov(bubb.f,x0,f_tol=5e-8,method="lgmres",verbose=True)
    #y = scipy.optimize.newton(bubb.f,x0) 
    y = bubb.solve(x0)
    #y = scipy.optimize.fsolve(bubb.f,x0)
    
    x = (np.dot(bubb.p1(),y) + np.dot(bubb.p2(),bubb.x_geg())).reshape(-1,4)
    X = np.append(X,x[init.ns+1::,:],axis=0)
    print(x[2*init.ns+1,3])
    #h = bubb.riemann_pos(x)

    if init.integration == "global":
        break
    # """ plot initial characteristic net  """
    # for j in range(0,len(init.elspos()[:,0])):
    #     mpos = init.elspos()[j,0]
    #     npos = init.elspos()[j,1]
    #     mneg = init.elsneg()[j,0]
    #     nneg = init.elsneg()[j,1]
    #     plt.plot([x[mpos,0],x[npos,0]],[x[mpos,1],x[npos,1]],"-",
    #                          color="tab:blue",linewidth=0.5)
    #     plt.plot([x[mneg,0],x[nneg,0]],[x[mneg,1],x[nneg,1]],"-",
    #                          color="tab:blue",linewidth=0.5)
    # plt.show()
    """ update """
    x = bubb.update(x)
    tf = x[3*init.ns+1,1]
    i += 1

    G = np.append(G,np.array([[x[3*init.ns+1,1],x[3*init.ns+1,2]]]),axis=0)
    # """ plot initial characteristic net  """
    # for j in range(0,len(init.elspos()[:,0])):
    #     mpos = init.elspos()[j,0]
    #     npos = init.elspos()[j,1]
    #     mneg = init.elsneg()[j,0]
    #     nneg = init.elsneg()[j,1]
    #     plt.plot([x[mpos,0],x[npos,0]],[x[mpos,1],x[npos,1]],"-",
    #                          color="tab:blue",linewidth=0.5)
    #     plt.plot([x[mneg,0],x[nneg,0]],[x[mneg,1],x[nneg,1]],"-",
    #                          color="tab:blue",linewidth=0.5)
    # plt.show()
G = np.delete(G,0,0)

np.savetxt("./data/data_error_analysis/q_s1.dat",G)
np.savetxt("./data/data_solution/X.dat",X)
np.savetxt("./data/data_solution/elspos.dat",init.elspos(),fmt='%i')
np.savetxt("./data/data_solution/elsneg.dat",init.elsneg(),fmt='%i')


""" post processing """
timeslabs = int((len(X)-(init.ns+1))/((3*init.ns+2)-(init.ns+1)))
qs0 = []
qs1 = []
as0 = []
as1 = []
for i in range(0,timeslabs):
    xi = X[i*(2*init.ns+1):i*(2*init.ns+1)+3*init.ns+2,:]
    p1 = 0
    p2 = 2*init.ns+1
    p3 = init.ns
    p4 = p3+2*init.ns+1
    qs0 += [xi[p1,1],xi[p1,2]]
    qs1 += [xi[p3,1],xi[p3,2]]
    as0 += [xi[p1,1],xi[p1,3]]
    as1 += [xi[p3,1],xi[p3,3]]
qs0 = np.array(qs0).reshape(-1,2)
qs1 = np.array(qs1).reshape(-1,2)
as0 = np.array(as0).reshape(-1,2)
as1 = np.array(as1).reshape(-1,2)
np.savetxt("./data/data_solution/qs0.dat",qs0)
np.savetxt("./data/data_solution/qs1.dat",qs1)
np.savetxt("./data/data_solution/as0.dat",as0)
np.savetxt("./data/data_solution/as1.dat",as1)
