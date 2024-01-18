import numpy as np
import matplotlib.pyplot as plt

data_xs, data_ys100, data_ks, data_taus, data_yas100 = np.loadtxt("./output_simple_100.dat", comments='#', unpack=True)
data_xc, data_yc100, data_kc, data_tauc, data_yac100= np.loadtxt("./output_cons_100.dat", comments='#', unpack=True)
data_xs, data_ys1, data_ks, data_taus, data_yas1 = np.loadtxt("./output_simple_1.dat", comments='#', unpack=True)
data_xc, data_yc1, data_kc, data_tauc, data_yac1= np.loadtxt("./output_cons_1.dat", comments='#', unpack=True)

ngrid = 1000

#print(data_x[4])
fig = plt.figure(figsize=(12, 12))
font_dict = dict(size=32)

ax1 = fig.add_subplot(221)
ax1.plot(data_xs[ngrid*1:ngrid*2], data_ys100[ngrid*1:ngrid*2], "-", lw=3,color="k", label="t=1.5")
ax1.plot(data_xs[ngrid*2:ngrid*3], data_ys100[ngrid*2:ngrid*3], "--", lw=3,color="k", label="t=3.0")
ax1.plot(data_xs[ngrid*1:ngrid*2], data_yas100[ngrid*1:ngrid*2], "-", lw=1,color="k", label="t=1.5(analytic)")
ax1.plot(data_xs[ngrid*2:ngrid*3], data_yas100[ngrid*2:ngrid*3], "--", lw=1,color="k", label="t=3.0(analytic)")
ax1.set_xlim(0, 100)
#ax.set_ylim(0, 10)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.legend(loc="lower right")
ax1.set_title(r"$\alpha = 100$, simple", size = 20)
ax1.text(70, 0.5, "(1)", fontdict=font_dict)

ax3 = fig.add_subplot(223)
ax3.plot(data_xc[ngrid*1:ngrid*2], data_yc100[ngrid*1:ngrid*2], "-", lw=3, color="k", label="t=1.5")
ax3.plot(data_xc[ngrid*2:ngrid*3], data_yc100[ngrid*2:ngrid*3], "--", lw=3, color="k", label="t=3.0")
ax3.plot(data_xc[ngrid*1:ngrid*2], data_yac100[ngrid*1:ngrid*2], "-", lw=1, color="k", label="t=1.5(analytic)")
ax3.plot(data_xc[ngrid*2:ngrid*3], data_yac100[ngrid*2:ngrid*3], "--", lw=1, color="k", label="t=3.0(analytic)")

ax3.set_xlim(0, 100)
ax3.set_xlabel("x")
ax3.set_ylabel("y")
ax3.legend(loc="lower right")
ax3.set_title(r"$\alpha = 100$, conserve", size = 20)
ax3.text(70, 0.5, "(3)", fontdict=font_dict)


ax2 = fig.add_subplot(222)
ax2.plot(data_xs[ngrid*1:ngrid*2], data_ys1[ngrid*1:ngrid*2], "-", lw=3, color="k", label="t=1.5")
ax2.plot(data_xs[ngrid*2:ngrid*3], data_ys1[ngrid*2:ngrid*3], "--", lw=3, color="k", label="t=3.0")
ax2.plot(data_xs[ngrid*1:ngrid*2], data_yas1[ngrid*1:ngrid*2], "-", lw=1, color="k", label="t=1.5(analytic)")
ax2.plot(data_xs[ngrid*2:ngrid*3], data_yas1[ngrid*2:ngrid*3], "--", lw=1, color="k", label="t=3.0(analytic)")

ax2.set_xlim(0, 100)
#ax.set_ylim(0, 10)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.legend(loc="lower right")
ax2.set_title(r"$\alpha = 1$, simple", size = 20)
ax2.text(70, 0.5, "(2)", fontdict=font_dict)


ax4 = fig.add_subplot(224)
ax4.plot(data_xc[ngrid*1:ngrid*2], data_yc1[ngrid*1:ngrid*2], "-", lw=3, color="k", label="t=1.5")
ax4.plot(data_xc[ngrid*2:ngrid*3], data_yc1[ngrid*2:ngrid*3], "--", lw=3, color="k", label="t=3.0")
ax4.plot(data_xc[ngrid*1:ngrid*2], data_yac1[ngrid*1:ngrid*2], "-", lw=1, color="k", label="t=1.5(analytic)")
ax4.plot(data_xc[ngrid*2:ngrid*3], data_yac1[ngrid*2:ngrid*3], "--", lw=1, color="k", label="t=3.0(analytic)")

ax4.set_xlim(0, 100)
#ax.set_ylim(0, 10)
ax4.set_xlabel("x")
ax4.set_ylabel("y")
ax4.legend(loc="lower right")
ax4.set_title(r"$\alpha = 1$, conserve", size = 20)
ax4.text(70, 0.5, "(4)", fontdict=font_dict)

plt.savefig("photon_cons.png")
plt.show()
