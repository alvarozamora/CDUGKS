import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import struct
import glob
import argparse
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=int, default=1)
parser.add_argument('-n', type=int, default=64)

args = parser.parse_args()


files = glob.glob('Data/rho_*')
vxfiles = glob.glob('Data/rhovx_*')
vyfiles = glob.glob('Data/rhovy_*')
vzfiles = glob.glob('Data/rhovz_*')
Efiles = glob.glob('Data/rhoE_*')

size = 8
num  = args.n
#num  = 32
type = 'd' #d is double, f is float, i is integer

testproblem = args.p

files.sort()
vxfiles.sort()
vyfiles.sort()
vzfiles.sort()
Efiles.sort()
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


if testproblem == 5:
	'''
	for file in vxfiles:
		n = file[11:15]
		n = file[11:]
		f = open(file, 'rb')

		X = f.read(num*num*size)
		X = np.array(struct.unpack(type*num*num, X))
		N = int(np.round(np.sqrt(len(X))))
		X = X.reshape((N,N))

		rhomin = np.minimum(X.min(),rhomin)
		rhomax = np.maximum(X.max(),rhomax)

	print(f'Max and Min density are {rhomax:.3f} and {rhomin:.3f}')
	'''

	f0 = open(vxfiles[0], 'rb')
	rhofile0 = open(files[0],'rb')
	X0 = f0.read(num*num*size)
	X0 = np.array(struct.unpack(type*num*num, X0))
	X0 = X0.reshape((args.n, args.n))
	rho0 = rhofile0.read(num*num*size)
	rho0 = np.array(struct.unpack(type*num*num, rho0))
	rho0 = rho0.reshape((args.n, args.n))
	#pdb.set_trace()
	V0 = X0/rho0
	V0 = V0[:,0]
	for q, file in enumerate(vxfiles,0):
		n = file[9:13]
		n = file[8:]
		f = open(file, 'rb')

		X = f.read(num*num*size)
		X = np.array(struct.unpack(type*num*num, X))

		X = X.reshape((args.n, args.n))

		rhofile = open(f'Data/rho_{q:04d}','rb')
		rho = rhofile.read(num*num*size)
		rho = np.array(struct.unpack(type*num*num, rho))
		rho = rho.reshape((args.n, args.n))

		X = X/rho
		#plt.imshow(X)
		#plt.savefig('Check/VXmap'+n+'.png')
		#print('Check/VXmap'+n+'.png')

		plt.cla()
		plt.clf()


		VX = X[:,0]
		y = np.linspace(1/2/len(VX), 1 - 1/2/len(VX), len(VX))
		plt.plot(y, V0,'k--',alpha=0.6, label='ICs')
		plt.plot(y, VX, label='Current')
		plt.ylim(-0.3, 0.3)
		plt.grid(alpha=0.2)
		plt.title(f'Horizontal Velocity (t = {q*4.0/401:.2f})')
		plt.legend()
		plt.savefig('Check/VXprof'+n+'.png')
		print('Check/VXprof'+n+'.png')

		plt.cla()
		plt.clf()
