import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('agg')
import matplotlib.pyplot as plt
import struct
import glob
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-p', type=int, default=1)
parser.add_argument('-n', type=int, default=64)
args = parser.parse_args()

files = glob.glob('Data/rho_*')
size = 8
num  = args.n

type = 'd' #d is double, f is float, i is integer

testproblem = args.p

files.sort()
files.reverse()
plt.figure()
if testproblem == 1:

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
		f.close()

if testproblem == 2:
	rhomin = 1.5
	rhomax = 1.5

	for file in files:
		n = file[9:13]
		n = file[8:]
		f = open(file, 'rb')

		X = f.read(num*num*size)
		X = np.array(struct.unpack(type*num*num, X))
		N = int(np.round(np.sqrt(len(X))))
		X = X.reshape((N,N))

		rhomin = np.minimum(X.min(),rhomin)
		rhomax = np.maximum(X.max(),rhomax)

	print(f'Max and Min density are {rhomax:.3f} and {rhomin:.3f}')
	#rhomin = 0.5
	#rhomax = 2.5

	for file in files:
		n = file[9:13]
		n = file[8:]
		f = open(file, 'rb')

		X = f.read(num*num*size)
		X = np.array(struct.unpack(type*num*num, X))
		N = int(np.round(np.sqrt(len(X))))
		X = X.reshape((N,N))

		grid = np.linspace(1/2/N, 1 - 1/2/N, N)
		xx, yy = np.meshgrid(grid, grid)
		plt.imshow(X, extent = (0,1,0,1))
		#plt.contourf(xx, yy, X, levels = np.linspace(rhomin,rhomax,300))
		plt.clim(0.5,2.5)
		plt.colorbar(ticks = [0.5, 1.0, 1.5, 2.0, 2.5])
		plt.savefig('Check/Rho'+n+'.png')
		print('Check/Rho'+n+'.png')

		plt.cla()
		plt.clf()

if testproblem == 3:
        for file in vxfiles:
                n = file[9:13]
                n = file[8:]
                f = open(file, 'rb')

                X = f.read(num*2*size)
                X = np.array(struct.unpack(type*num*2, X))

                X = X.reshape((args.n, 2))

                plt.imshow(X)
                plt.savefig('Check/VXmap'+n+'.png')
                print('Check/VXmap'+n+'.png')

                plt.cla()
                plt.clf()

                VX = X[:,0]
                y = np.linspace(1/2/len(VX), 1 - 1/2/len(VX), len(VX))
                plt.plot(y, VX)
                plt.ylim(-0.5, 0.5)
                plt.grid(alpha=0.2)
                plt.savefig('Check/VXprof'+n+'.png')
                print('Check/VXprof'+n+'.png')

                plt.cla()
                plt.clf()
