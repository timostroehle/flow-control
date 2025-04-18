import numpy as np
import matplotlib.pyplot as plt

q2_mid = np.loadtxt("./data/data_error_analysis/q_s1_2_mid.txt")
q4_mid = np.loadtxt("./data/data_error_analysis/q_s1_4_mid.txt")
q6_mid = np.loadtxt("./data/data_error_analysis/q_s1_6_mid.txt")
q8_mid = np.loadtxt("./data/data_error_analysis/q_s1_8_mid.txt")
q10_mid = np.loadtxt("./data/data_error_analysis/q_s1_10_mid.txt")
q12_mid = np.loadtxt("./data/data_error_analysis/q_s1_12_mid.txt")

q2_imp = np.loadtxt("./data/data_error_analysis/q_s1_2_imp.txt")
q4_imp = np.loadtxt("./data/data_error_analysis/q_s1_4_imp.txt")
q6_imp = np.loadtxt("./data/data_error_analysis/q_s1_6_imp.txt")
q8_imp = np.loadtxt("./data/data_error_analysis/q_s1_8_imp.txt")
q10_imp = np.loadtxt("./data/data_error_analysis/q_s1_10_imp.txt")
q12_imp = np.loadtxt("./data/data_error_analysis/q_s1_12_imp.txt")

q2_at_10_mid = np.interp([10],q2_mid[:,0],q2_mid[:,1])[0]
q4_at_10_mid = np.interp([10],q4_mid[:,0],q4_mid[:,1])[0]
q6_at_10_mid = np.interp([10],q6_mid[:,0],q6_mid[:,1])[0]
q8_at_10_mid = np.interp([10],q8_mid[:,0],q8_mid[:,1])[0]
q10_at_10_mid = np.interp([10],q10_mid[:,0],q10_mid[:,1])[0]
q12_at_10_mid = np.interp([10],q12_mid[:,0],q12_mid[:,1])[0]

q2_at_10_imp = np.interp([10],q2_imp[:,0],q2_imp[:,1])[0]
q4_at_10_imp = np.interp([10],q4_imp[:,0],q4_imp[:,1])[0]
q6_at_10_imp = np.interp([10],q6_imp[:,0],q6_imp[:,1])[0]
q8_at_10_imp = np.interp([10],q8_imp[:,0],q8_imp[:,1])[0]
q10_at_10_imp = np.interp([10],q10_imp[:,0],q10_imp[:,1])[0]
q12_at_10_imp = np.interp([10],q12_imp[:,0],q12_imp[:,1])[0]

plt.plot(q2_mid[:,0],q2_mid[:,1],"-",color="magenta")
plt.plot(10,q2_at_10_mid,"o")
plt.plot(q4_mid[:,0],q4_mid[:,1],"-",color="red")
plt.plot(10,q4_at_10_mid,"o")
plt.plot(q6_mid[:,0],q6_mid[:,1],"-",color="blue")
plt.plot(10,q6_at_10_mid,"o")
plt.plot(q8_mid[:,0],q8_mid[:,1],"-",color="green")
plt.plot(10,q8_at_10_mid,"o")
plt.plot(q10_mid[:,0],q10_mid[:,1],"-",color="orange")
plt.plot(10,q10_at_10_mid,"o")
plt.plot(q12_mid[:,0],q12_mid[:,1],"-",color="orange")
plt.plot(10,q12_at_10_mid,"o")
plt.show()

plt.plot(q2_imp[:,0],q2_imp[:,1],"-",color="magenta")
plt.plot(10,q2_at_10_imp,"o")
plt.plot(q4_imp[:,0],q4_imp[:,1],"-",color="red")
plt.plot(10,q4_at_10_imp,"o")
plt.plot(q6_imp[:,0],q6_imp[:,1],"-",color="blue")
plt.plot(10,q6_at_10_imp,"o")
plt.plot(q8_imp[:,0],q8_imp[:,1],"-",color="green")
plt.plot(10,q8_at_10_imp,"o")
plt.plot(q10_imp[:,0],q10_imp[:,1],"-",color="green")
plt.plot(10,q10_at_10_imp,"o")
plt.plot(q12_imp[:,0],q12_imp[:,1],"-",color="green")
plt.plot(10,q12_at_10_imp,"o")
plt.show()

q_ana = 0.0107584 # 0.010687187983212611

err2_mid = abs(q2_at_10_mid-q_ana)/abs(q_ana)
err4_mid = abs(q4_at_10_mid-q_ana)/abs(q_ana)
err6_mid = abs(q6_at_10_mid-q_ana)/abs(q_ana)
err8_mid = abs(q8_at_10_mid-q_ana)/abs(q_ana)
err10_mid = abs(q10_at_10_mid-q_ana)/abs(q_ana)
err12_mid = abs(q12_at_10_mid-q_ana)/abs(q_ana)
err_mid = np.array([[2,err2_mid],
                    [4,err4_mid],
                    [6,err6_mid],
                    [8,err8_mid],
                    [10,err10_mid],
                    [12,err12_mid]])

err2_imp = abs(q2_at_10_imp-q_ana)/abs(q_ana)
err4_imp = abs(q4_at_10_imp-q_ana)/abs(q_ana)
err6_imp = abs(q6_at_10_imp-q_ana)/abs(q_ana)
err8_imp = abs(q8_at_10_imp-q_ana)/abs(q_ana)
err10_imp = abs(q10_at_10_imp-q_ana)/abs(q_ana)
err12_imp = abs(q12_at_10_imp-q_ana)/abs(q_ana)
err_imp = np.array([[2,err2_imp],
                    [4,err4_imp],
                    [6,err6_imp],
                    [8,err8_imp],
                    [10,err10_imp],
                    [12,err12_imp]])

E_mid = np.append(np.array([[2,4,6,8,10,12]]),
                  np.array([[err2_mid,err4_mid,err6_mid,err8_mid,err10_mid,err12_mid]]),
                  axis=0).T
E_imp = np.append(np.array([[2,4,6,8,10,12]]),
                  np.array([[err2_imp,err4_imp,err6_imp,err8_imp,err10_imp,err12_imp]]),
                  axis=0).T

np.savetxt("./data/data_error_analysis/err_mid.dat",E_mid)
np.savetxt("./data/data_error_analysis/err_imp.dat",E_imp)

plt.loglog([2,4,6,8,10,12],[err2_mid,err4_mid,err6_mid,err8_mid,err10_mid,err12_mid],"-",
           color="tab:blue",label="Mid-point")
plt.loglog([2,4,6,8,10,12],[err2_mid,err4_mid,err6_mid,err8_mid,err10_mid,err12_mid],"o",
           color="tab:blue")
plt.loglog([2e0,4e0,2e0,2e0],[1e-5,1e-5,4*1e-5,1e-5],"-",color="tab:blue")
plt.text(2.1e0,1.7e-5,"2")

plt.loglog([2,4,6,8,10,12],[err2_imp,err4_imp,err6_imp,err8_imp,err10_imp,err12_imp],"-",
           color="tab:orange",label="Euler impl.")
plt.loglog([2,4,6,8,10,12],[err2_imp,err4_imp,err6_imp,err8_imp,err10_imp,err12_imp],"o",
           color="tab:orange")
plt.loglog([2e0,4e0,2e0,2e0],[1e-2,1e-2,2*1e-2,1e-2],"-",color="tab:orange")
plt.text(2.1e0,1.2e-2,"1")

plt.grid(True, which="both", ls="-")
plt.ylabel("\( \\varepsilon(n_{s}) \)")
plt.xlabel("\( n_{s} \)")
plt.legend(loc="center right")
plt.savefig('moc_num_error_analysis.pgf', format='pgf',bbox_inches="tight")
plt.show()
