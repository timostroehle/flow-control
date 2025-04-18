import numpy as np
import matplotlib.pyplot as plt

q2_mid = np.loadtxt("./q_s1_mp2.dat")
q4_mid = np.loadtxt("./q_s1_mp4.dat")
q8_mid = np.loadtxt("./q_s1_mp8.dat")
q16_mid = np.loadtxt("./q_s1_mp16.dat")
q24_mid = np.loadtxt("./q_s1_mp24.dat")
q32_mid = np.loadtxt("./q_s1_mp32.dat")

q2_imp = np.loadtxt("./q_s1_ie2.dat")
q4_imp = np.loadtxt("./q_s1_ie4.dat")
q8_imp = np.loadtxt("./q_s1_ie8.dat")
q16_imp = np.loadtxt("./q_s1_ie16.dat")
q24_imp = np.loadtxt("./q_s1_ie24.dat")
q32_imp = np.loadtxt("./q_s1_ie32.dat")

p = 0.08
q2_at_p_mid = np.interp([p],q2_mid[:,0],q2_mid[:,1])[0]
q4_at_p_mid = np.interp([p],q4_mid[:,0],q4_mid[:,1])[0]
q8_at_p_mid = np.interp([p],q8_mid[:,0],q8_mid[:,1])[0]
q16_at_p_mid = np.interp([p],q16_mid[:,0],q16_mid[:,1])[0]
q24_at_p_mid = np.interp([p],q24_mid[:,0],q24_mid[:,1])[0]
q32_at_p_mid = np.interp([p],q32_mid[:,0],q32_mid[:,1])[0]

q2_at_p_imp = np.interp([p],q2_imp[:,0],q2_imp[:,1])[0]
q4_at_p_imp = np.interp([p],q4_imp[:,0],q4_imp[:,1])[0]
q8_at_p_imp = np.interp([p],q8_imp[:,0],q8_imp[:,1])[0]
q16_at_p_imp = np.interp([p],q16_imp[:,0],q16_imp[:,1])[0]
q24_at_p_imp = np.interp([p],q24_imp[:,0],q24_imp[:,1])[0]
q32_at_p_imp = np.interp([p],q32_imp[:,0],q32_imp[:,1])[0]


plt.plot(q2_mid[:,0],q2_mid[:,1],"-",color="magenta")
plt.plot(p,q2_at_p_mid,"o")
plt.plot(q4_mid[:,0],q4_mid[:,1],"-",color="red")
plt.plot(p,q4_at_p_mid,"o")
plt.plot(q8_mid[:,0],q8_mid[:,1],"-",color="green")
plt.plot(p,q8_at_p_mid,"o")
plt.plot(q16_mid[:,0],q16_mid[:,1],"-",color="orange")
plt.plot(p,q16_at_p_mid,"o")
plt.plot(q24_mid[:,0],q24_mid[:,1],"-",color="orange")
plt.plot(p,q24_at_p_mid,"o")
plt.plot(q32_mid[:,0],q32_mid[:,1],"-",color="orange")
plt.plot(p,q32_at_p_mid,"o")
plt.show()

plt.plot(q2_imp[:,0],q2_imp[:,1],"-",color="magenta")
plt.plot(p,q2_at_p_imp,"o")
plt.plot(q4_imp[:,0],q4_imp[:,1],"-",color="red")
plt.plot(p,q4_at_p_imp,"o")
plt.plot(q8_imp[:,0],q8_imp[:,1],"-",color="green")
plt.plot(p,q8_at_p_imp,"o")
plt.plot(q16_imp[:,0],q16_imp[:,1],"-",color="green")
plt.plot(p,q16_at_p_imp,"o")
plt.plot(q24_imp[:,0],q24_imp[:,1],"-",color="green")
plt.plot(p,q24_at_p_imp,"o")
plt.plot(q32_imp[:,0],q32_imp[:,1],"-",color="green")
plt.plot(p,q32_at_p_imp,"o")
plt.show()

q_ana = 1.701398691981323708e-01 #1.699995597788319657e-01 #0.0107584 #0.010687187983212611

print(q_ana)
print(q32_at_p_mid)

err2_mid = abs(q2_at_p_mid-q_ana)/abs(q_ana)
err4_mid = abs(q4_at_p_mid-q_ana)/abs(q_ana)
err8_mid = abs(q8_at_p_mid-q_ana)/abs(q_ana)
err16_mid = abs(q16_at_p_mid-q_ana)/abs(q_ana)
err24_mid = abs(q24_at_p_mid-q_ana)/abs(q_ana)
err32_mid = abs(q32_at_p_mid-q_ana)/abs(q_ana)
err_mid = np.array([[2,err2_mid],
                    [4,err4_mid],
                    [8,err8_mid],
                    [16,err16_mid],
                    [24,err24_mid],
                    [32,err32_mid]])

err2_imp = abs(q2_at_p_imp-q_ana)/abs(q_ana)
err4_imp = abs(q4_at_p_imp-q_ana)/abs(q_ana)
err8_imp = abs(q8_at_p_imp-q_ana)/abs(q_ana)
err16_imp = abs(q16_at_p_imp-q_ana)/abs(q_ana)
err24_imp = abs(q24_at_p_imp-q_ana)/abs(q_ana)
err32_imp = abs(q32_at_p_imp-q_ana)/abs(q_ana)
err_imp = np.array([[2,err2_imp],
                    [4,err4_imp],
                    [8,err8_imp],
                    [16,err16_imp],
                    [24,err24_imp],
                    [32,err32_imp]])

E_mid = np.append(np.array([[2,4,8,16,24,32]]),
                  np.array([[err2_mid,err4_mid,err8_mid,err16_mid,err24_mid,err32_mid]]),
                  axis=0).T
E_imp = np.append(np.array([[2,4,8,16,24,32]]),
                  np.array([[err2_imp,err4_imp,err8_imp,err16_imp,err24_imp,err32_imp]]),
                  axis=0).T

np.savetxt("./err_mid.dat",E_mid)
np.savetxt("./err_imp.dat",E_imp)

plt.loglog([2,4,8,16,24,32],[err2_mid,err4_mid,err8_mid,err16_mid,err24_mid,err32_mid],"-",
           color="tab:blue",label="Mid-point")
plt.loglog([2,4,8,16,24,32],[err2_mid,err4_mid,err8_mid,err16_mid,err24_mid,err32_mid],"o",
           color="tab:blue")
plt.loglog([2e0,4e0,2e0,2e0],[1e-6,1e-6,4*1e-6,1e-6],"-",color="tab:blue")
plt.text(2.1e0,1.7e-6,"2")

plt.loglog([2,4,8,16,24,32],[err2_imp,err4_imp,err8_imp,err16_imp,err24_imp,err32_imp],"-",
           color="tab:orange",label="Euler impl.")
plt.loglog([2,4,8,16,24,32],[err2_imp,err4_imp,err8_imp,err16_imp,err24_imp,err32_imp],"o",
           color="tab:orange")
plt.loglog([2e0,4e0,2e0,2e0],[1e-2,1e-2,2*1e-2,1e-2],"-",color="tab:orange")
plt.text(2.1e0,1.2e-2,"1")

plt.grid(True, which="both", ls="-")
plt.ylabel("\( \\varepsilon(n_{s}) \)")
plt.xlabel("\( n_{s} \)")
plt.legend(loc="center right")
plt.savefig('moc_num_error_analysis.pgf', format='pgf',bbox_inches="tight")
plt.show()
