import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import struct
import glob

files = glob.glob('Data/rho*')
size = 16
num  = 256
type = 'd' #d is double, f is float, i is integer

files.sort()
plt.figure()
for file in files:
	n = file[9:13]
	f = open(file, 'rb')

	X = f.read(num*size)
	X = np.array(struct.unpack(type*num, X))

	plt.xlim(0,1)
	plt.ylim(0,1.2)
	plt.plot(np.linspace(0.5/len(X),1-0.5/len(X), len(X)), X, 'o-')
	plt.grid()
	plt.savefig('Check/Rho'+n+'.png')
	print('Check/Rho'+n+'.png')

	plt.cla()
	plt.clf()
