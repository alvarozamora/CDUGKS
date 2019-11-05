import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import struct
import glob

files = glob.glob('Data2/rho*')
size = 8
num  = 256
type = 'd' #d is double, f is float, i is integer

testproblem = 1

files.sort()
plt.figure()
if testproblem == 1:

	for file in files:
		n = file[10:14]
		f = open(file, 'rb')

		X = f.read(num*size)
		X = np.array(struct.unpack(type*num, X))

		plt.xlim(0,1)
		plt.ylim(0,1.2)
		plt.plot(np.linspace(0.5/len(X),1-0.5/len(X), len(X)), X, 'o-')
		plt.grid()
		plt.savefig('Check2/Rho'+n+'.png')
		print('Check2/Rho'+n+'.png')

		plt.cla()
		plt.clf()
		f.close()

if testproblem == 2:

	for file in files:
		n = file[9:13]
		f = open(file, 'rb')

		X = f.read(num*num*size)
		X = np.array(struct.unpack(type*num*num, X))
		N = int(np.round(np.sqrt(len(X))))
		X = X.reshape((N,N))

		plt.imshow(X)
		plt.savefig('Check/Rho'+n+'.png')
		print('Check/Rho'+n+'.png')

		plt.cla()
		plt.clf()
