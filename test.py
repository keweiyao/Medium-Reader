import medium
import matplotlib.pyplot as plt
import sys
import h5py
import numpy as np

f = h5py.File(sys.argv[1], 'r')
A = medium.Medium(sys.argv[1])
t = []
NQ = 10
info = np.zeros([15, NQ, 300])

pos = (2.*np.random.rand(NQ*2).reshape(NQ,2) - 1.)*4.
keys = ['Temp', 'Vx', 'Vy', 'e', 'p', 'pi00','pi01','pi02','pi03','pi11','pi12','pi13','pi22','pi23','pi33']
labels = [r'$T$ [GeV]', r'$v_x$', r'$v_y$', r'$e$ [GeV/fm${}^3$]', r'$p$ [GeV/fm${}^3$]', r'$\pi^{00}/(e+p)$',r'$\pi^{01}/(e+p)$',r'$\pi^{02}/(e+p)$',r'$\pi^{03}/(e+p)$',r'$\pi^{11}/(e+p)$',r'$\pi^{12}/(e+p)$',r'$\pi^{13}/(e+p)$',r'$\pi^{22}/(e+p)$',r'$\pi^{23}/(e+p)$',r'$\pi^{33}/(e+p)$']
for i in range(300):
	t.append(0.6+i*0.04+0.0)
	A.load_next()
	for j, x in enumerate(pos):
		for k, key in enumerate(keys):
			info[k, j, i] = A.interpF(t[i], x, keys[k])
info[5:] /= (info[3] + info[4])

plt.figure(figsize=(15,15))
for i, data in enumerate(info):
	plt.subplot(5,3,i+1)
	plt.ylabel(labels[i], size=20)
	for q in data:
		plt.plot(t, q, color = 'r', alpha = 0.3)
	plt.subplots_adjust(wspace=0.3, hspace=0.25)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.xlabel(r'$\tau$ [fm/c]' if i>11 else '', size=20)
	plt.xticks([0, 5, 10, 15] if i>11 else [])

plt.show()

