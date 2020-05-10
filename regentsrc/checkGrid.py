import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('agg')
import matplotlib.pyplot as plt
import struct
import glob
import argparse
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=int, default=0)
parser.add_argument('-n', type=int, default=64)
args = parser.parse_args()

rhofiles = glob.glob('Data/rho_*')
vxfiles = glob.glob('Data/rhovx_*')
vyfiles = glob.glob('Data/rhovy_*')
vzfiles = glob.glob('Data/rhovz_*')
Efiles = glob.glob('Data/rhoE_*')

phasefiles = glob.glob('Data/phase_*')


size = 8
num  = args.n

type = 'd' #d is double, f is float, i is integer

testproblem = args.p

rhofiles.sort()
rhofiles.reverse()
vxfiles.sort()
vxfiles.reverse()
vyfiles.sort()
vyfiles.reverse()
vzfiles.sort()
vzfiles.reverse()
Efiles.sort()
Efiles.reverse()
phasefiles.sort()
#phasefiles.reverse()

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

if testproblem == 6:
	from matplotlib.ticker import NullFormatter
	for file in phasefiles:
		pass
		'''
		n = file[-5:]
		f = open(file, 'rb')

		F = f.read(num*num*size)
		F = np.array(struct.unpack(type*num*num,F))


		F = F.reshape((num,num))
		plt.imshow(F,origin='lower')
		plt.savefig("Check/phase"+n+".png")
		plt.figure()
		avg1 = F.mean(axis=0)
		avg2 = F.mean(axis=1)

		v = np.linspace(-10,10,num)
		plt.plot(v, avg2)
		plt.savefig('Check/phaseavg'+n+'.png')
		print(f"Problem 6: Finished {n[1:]}")
		'''
	for file in phasefiles:
		n = file[-5:]
		f = open(file, 'rb')

		F = f.read(num*num*size)
		F = np.array(struct.unpack(type*num*num,F))

		x = np.linspace(0, 1,num)
		v = np.linspace(-10,10,num)

		F = F.reshape((num,num))


		nullfmt = NullFormatter()         # no labels

		# definitions for the axes
		left, width = 0.15, 0.6
		bottom, height = 0.15, 0.6
		bottom_h = left_h = left + width + 0.05

		rect_Phase = [left, bottom, width, height]
		rect_avgx = [left, bottom_h, width, 0.15]
		rect_avgy = [left_h, bottom, 0.15, height]

		# start with a rectangular Figure
		plt.figure(1, figsize=(8, 8))
		plt.clf()

		print(rect_Phase)
		axPhase = plt.axes(rect_Phase)
		axavgx = plt.axes(rect_avgx)
		axavgy = plt.axes(rect_avgy)

		# no labels
		axavgx.xaxis.set_major_formatter(nullfmt)
		axavgy.yaxis.set_major_formatter(nullfmt)

		# the Phase plot:
		axPhase.imshow(F, extent=[0,1,-10,10],aspect='auto')#,aspect=1/20)
		axPhase.set_ylabel('Velocity')
		axPhase.set_xlabel('Position')
		axPhase.set_ylim((-10, 10))
		axPhase.set_xlim((0,1))

		#pdb.set_trace()
		axavgx.plot(x, F.mean(axis=0))
		axavgy.plot(F.mean(axis=1), v)#, orientation='horizontal')

		axavgx.set_xlim(axPhase.get_xlim())
		#axavgx.set_ylim(0,0.1)

		#axavgy.set_xlim(0,1.5)
		axavgy.set_ylim(axPhase.get_ylim())

		plt.savefig("Check/phase"+n+'.png')
		print(f"Problem 6: Finished {n[1:]}")
