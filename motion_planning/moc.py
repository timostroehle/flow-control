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
        #self.da = 0.05
        self.dq = 2.0
        self.method = "mp" # for midpoint-rule or "ie" for implicit euler
        self.problem = "inv"  #"dir" or  "inv"
        """ loading data """
        self.path_to_data_init = "./data/data_init/"
        self.path_to_data_solution = "./data/data_solution/"
        self.path_to_data_error_analysis = "./data/data_error_analysis/"
        
        self.ns = int(np.loadtxt(self.path_to_data_init+"param.dat")[0])
        self.nt = int(np.loadtxt(self.path_to_data_init+"param.dat")[1])
        self.c0 = np.loadtxt(self.path_to_data_init+"param.dat")[2]
        self.rho= np.loadtxt(self.path_to_data_init+"param.dat")[3]
        self.a0 = np.loadtxt(self.path_to_data_init+"a0.dat")
        self.q0 = np.loadtxt(self.path_to_data_init+"q0.dat")
        self.elspos = np.loadtxt(self.path_to_data_init+"elspos.dat", dtype="int")
        self.elsneg = np.loadtxt(self.path_to_data_init+"elsneg.dat", dtype="int")
        self.x_init = np.loadtxt(self.path_to_data_init+"xinit.dat")
        self.nkn = len(self.x_init)
        self.bc = np.loadtxt(self.path_to_data_init+"bc.dat",dtype="int").tolist()
        self.bc_trj = np.loadtxt(self.path_to_data_init+"bc_trj.dat",dtype="int",ndmin=1).tolist()
    
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

    def update(self,x):
        xn = np.array(x)
        for i in range(0,2*init.ns+1):
            """ update """
            xn[i,0] = x[i+(2*init.ns+1),0]
            xn[i,1] = x[i+(2*init.ns+1),1]
            xn[i,2] = x[i+(2*init.ns+1),2]
            xn[i,3] = x[i+(2*init.ns+1),3]
            """ extrapolation """
            dti = x[i+(2*init.ns+1),1] - x[i,1]
            xn[i+(2*init.ns+1),0] = x[i+(2*init.ns+1),0]
            xn[i+(2*init.ns+1),1] = x[i+(2*init.ns+1),1] + dti 
            xn[i+(2*init.ns+1),2] = x[i+(2*init.ns+1),2]
            xn[i+(2*init.ns+1),3] = x[i+(2*init.ns+1),3]
        return xn
    
    def a_pre(self,t):
        t0 = 1.0; t1 = 4.0; dt = t1-t0; da = init.da
        if t<t0:
            a = init.a0
        elif t>t1:
            a = init.a0+da
        else:
            """ C2-function (cos(cos())) """
            f = 0.5* (np.cos((t-t0)*np.pi/dt + np.pi) + 1)
            a = da*(0.5*np.cos(f*np.pi+np.pi) + 0.5) + 1
            """ C0 linear """
            #a = (1/(tf-t0))*t - t0/(tf-t0)
            """ C-1 constant """
            #a = init.a0 + da
        return a

    def q_pre(self,t):
        t0 = 0.005; t1 = 0.025; dt = t1-t0; dq = init.dq
        if t<t0:
            q = init.q0
        elif t>t1:
            q = init.q0+dq
        else:
            """ C2-function (cos(cos())) """
            f = 0.5* (np.cos((t-t0)*np.pi/dt + np.pi) + 1)
            q = dq*(0.5*np.cos(f*np.pi+np.pi) + 0.5) 
            """ C0 linear """
            #q = dq*((1/(t1-t0))*t - t0/(t1-t0))
            """ C-1 constant """
            #q = init.q0 + dq
        return q

    def dq_pre(self,t):
        t0 = 0.005; t1 = 0.025; dt = t1-t0; a = init.dq
        if t<t0:
            dq = 0
        elif t>t1:
            dq = 0
        else:
            """ C2-function (cos(cos())) """
            f = 0.5* (np.cos((t-t0)*np.pi/dt + np.pi) + 1)
            df = -0.5*np.sin((t-t0)*np.pi/dt + np.pi) * np.pi/dt
            dq = -0.5*a*np.sin(f*np.pi + np.pi)*np.pi*df 
            """ C0 linear """
            #dq = a*(1/(t1-t0))
            """ C-1 constant """
            #dq = 0
        return dq
    
    def beta(self,s):
        """ gaussian distribution, maximum at s=0 """
        sig = 0.1
        s0 = 0.0
        b1 = 0.05*(1/10)*(1/sig)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) +0.2#+ 5
        """ gaussian distribution of stiffness """
        sig = 0.25
        s0 = 0.5
        b = -0.5*(1/10)*(1/sig)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) + 0.3 + b1  # 2.5
        """ constant distribution of stiffness """
        #b = 2
        return (1/(0+1))*b+ 0

    def d_beta(self,s):
        """ gaussian distribution, maximum at s=0 """
        sig = 0.1
        s0 = 0.0
        a1 = 0.05*(1/10) 
        db1 = -a1*(2/sig**3)*(s-s0)*np.exp(-1.0*(((s-s0)**2)/(sig**2))) 
        """ gaussian distribution of stiffness """
        sig = 0.25
        s0 = 0.5
        a = -0.5*(1/10)
        db = -a*(2/sig**3)*(s-s0)*np.exp(-1.0*(((s-s0)**2)/(sig**2)))  + db1
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
                f += [(qm-qn)/(tm-tn) - (qm/am - c)*(am-an)/(tm-tn) + K*(qm/am) + a_m*db*g*(1/init.rho)]
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
                f += [(qm-qn)/(tm-tn)-(q_mid/a_mid-c)*(am-an)/(tm-tn) + K*(q_mid/a_mid) + a_mid*db*g*(1/init.rho)]
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
                f += [(sm-sn)/(tm-tn) - (qm/am - c)]
                #kom
                f += [(qm-qn)/(tm-tn) - (qm/am + c)*(am-an)/(tm-tn) + K*(qm/am) + a_m*db*g*(1/init.rho)]
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
                f += [(qm-qn)/(tm-tn)-(q_mid/a_mid+c)*(am-an)/(tm-tn) + K*(q_mid/a_mid) + a_mid*db*g*(1/init.rho)]
            else:
                print("numerical method is not implemented")
            """ ........... """
        """ indirect imposition of dirichlet bc for a on s=0,t """
        for i in range(0,len(init.bc_trj)):
            if init.problem=="dir":
                f += [x[i*(2*init.ns+1),3] - self.a_pre(x[i*(2*init.ns+1),1])]
            elif init.problem=="inv":
                f += [x[init.bc_trj[i],2] - self.q_pre(x[init.bc_trj[i],1])]
                continue
            
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
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db )
                dr1[0,4*m+3] =(-(1/(tm-tn))*( (q_mid/a_mid) + c )
                               +(0.5*(q_mid/a_mid**2) - (init.rho**(-0.5))*0.125*(2**(-0.5))*(b**0.5)*(a_mid**(-0.75)))*(am-an)/(tm-tn)
                               - 0.5*K*(q_mid/a_mid**2)
                               + (1/init.rho)*0.25*a_mid**(0.5)*db + (1/init.rho)*0.5*g*db )
                dr = np.append(dr,dr1,axis=0)
            else:
                print("numerical method is not implemented")
        """ algebraic constraints """
        dr1 = np.zeros([1,4*init.nkn])
        dr1[0,-2] = 1.0
        dr1[0,-3] = -self.dq_pre(x[-1,1])
        dr = np.append(dr,dr1,axis=0)
        """ delete dirichlet bc """
        dr = np.delete(dr,0,0)
        dr = np.delete(dr,np.s_[init.bc],1)
        return dr

    def num_tan(self,res,x):
        h_num =  1e-12
        dr = np.zeros([len(x),len(x)])
        for jj in range(0,len(x)):
            xh = np.array(x)
            xh[jj,0] = x[jj,0] + h_num
            res_h = self.f(xh)
            dr[:,jj] = (res_h[:,0]-res[:,0])/h_num
        return dr

    def solve(self,x):
        res_norm = 1
        while res_norm > 1e-7:
            print(res_norm)
            """ ..............................
            history = open('history.txt', 'a')
            history.write(str(res_norm))
            history.write('\n')
            history.close()
            ..............................."""
            """ residual vector ..........."""
            res = self.f(x)
            """ residual jacobian ........."""
            dres = self.ana_tan(x)
            #dres = self.num_tan(res,x)
            #dres = self.sym_tan()
            """ residual norm ............."""
            res_norm = np.linalg.norm(res)
            """ solve ....................."""
            dx = np.linalg.solve(dres,res)
            """ update ...................."""
            x = x - dx
        return x
    
""" initialisation ..................................................................."""
init = Init()
x = init.x_init
X = np.zeros([1,4])
bubb = Bubbles(x)

As0 = np.array([[x[0,1],x[0,3]]])
As1 = np.array([[x[2*init.ns,1],x[2*init.ns,3]]])
Qs0 = np.array([[x[0,1],x[0,2]]])
Qs1 = np.array([[x[2*init.ns,1],x[2*init.ns,2]]])

D = np.array([[x[0,1],x[2*init.ns,1]-x[0,1]]])
ht = np.array([[0,0]])

n = init.nt-init.ns
for i in range(0,n):
    print("time-step:"+ str(i)+"/"+str(n))
    bubb = Bubbles(x)
    x0 = bubb.x_ges()
    """ Solve """
    y = bubb.solve(x0)
    
    x = (np.dot(bubb.p1(),y) + np.dot(bubb.p2(),bubb.x_geg())).reshape(-1,4)
    X = np.append(X,x,axis=0)

    for j in range(0,len(init.elspos[:,0])):
        mpos = init.elspos[j,0]
        npos = init.elspos[j,1]
        mneg = init.elsneg[j,0]
        nneg = init.elsneg[j,1]
        plt.plot([x[mpos,0],x[npos,0]],[x[mpos,1],x[npos,1]],"-",
                 color="tab:blue",linewidth=0.5)
        plt.plot([x[mneg,0],x[nneg,0]],[x[mneg,1],x[nneg,1]],"-",
                 color="tab:blue",linewidth=0.5)

    As0 = np.append(As0,np.array([[x[2*init.ns+1,1],x[2*init.ns+1,3]]]),axis=0)
    As1 = np.append(As1,np.array([[x[4*init.ns+1,1],x[4*init.ns+1,3]]]),axis=0)
    Qs0 = np.append(Qs0,np.array([[x[2*init.ns+1,1],x[2*init.ns+1,2]]]),axis=0)
    Qs1 = np.append(Qs1,np.array([[x[4*init.ns+1,1],x[4*init.ns+1,2]]]),axis=0)

    D = np.append(D,np.array([[x[0,1],x[2*init.ns,1]-x[0,1]]]),axis=0)
    ht = np.append(ht,np.array([[x[0,1],x[2*init.ns+1,1]-x[0,1]]]),axis=0)

    
    """ update """
    x = bubb.update(x)

X = np.delete(X,0,0)
ht = np.delete(ht,0,0)
np.savetxt(init.path_to_data_solution+"X.dat",X)
np.savetxt(init.path_to_data_solution+"as0.dat",As0)
np.savetxt(init.path_to_data_solution+"as1.dat",As1)
np.savetxt(init.path_to_data_solution+"qs0.dat",Qs0)
np.savetxt(init.path_to_data_solution+"qs1.dat",Qs1)
np.savetxt(init.path_to_data_solution+"d.dat",D)
np.savetxt(init.path_to_data_solution+"ht.dat",ht)
