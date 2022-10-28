import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('agg')
import matplotlib.pyplot as plt
import struct
import glob
import argparse
import pdb
import os
import io
import imageio as imo
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=int, default=0)
parser.add_argument('-n', type=int, default=64)
parser.add_argument('-v', type=int, default=0) # set below
args = parser.parse_args()
args.v = args.n if args.v == 0 else args.v

rhofiles = glob.glob('Data/rho_*')
vxfiles = glob.glob('Data/rhovx_*')
vyfiles = glob.glob('Data/rhovy_*')
vzfiles = glob.glob('Data/rhovz_*')
Efiles = glob.glob('Data/rhoE_*')

phasefiles = glob.glob('Data/phase_*')


size = 8
num  = args.n
numv = args.v

type = 'd' #d is double, f is float, i is integer

testproblem = args.p

rhofiles.sort()
#rhofiles.reverse()
vxfiles.sort()
#vxfiles.reverse()
vyfiles.sort()
vyfiles.reverse()
vzfiles.sort()
vzfiles.reverse()
Efiles.sort()
Efiles.reverse()
phasefiles.sort()
phasefiles.reverse()

os.makedirs("Plots", exist_ok=True)

plt.figure()
if testproblem == 1:

	for file in rhofiles[:]:
		n = file[9:13]
		f = open(file, 'rb')

		X = f.read(num*size)
		X = np.array(struct.unpack(type*num, X))

		plt.xlim(0,1)
		plt.ylim(0,1.2)
		plt.plot(np.linspace(0.5/len(X),1-0.5/len(X), len(X)), X)#, 'o-')
		plt.grid()
		plt.savefig('Plots/Rho'+n+'.png')
		print('Plots/Rho'+n+'.png')

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
		plt.imshow(X, origin='lower',extent = (0,1,0,1))
		#plt.contourf(xx, yy, X, levels = np.linspace(rhomin,rhomax,300))
		plt.clim(0.5,2.5)
		plt.colorbar(ticks = [0.5, 1.0, 1.5, 2.0, 2.5])
		plt.savefig('Plots/Rho'+n+'.png')
		print('Plots/Rho'+n+'.png')

		plt.cla()
		plt.clf()

#if testproblem == 3:
if testproblem == 5:
	#visc = np.array([1e-2, 1e-4, 1e-6])*(10/7)**0.81
	visc = 1e-2
	times = np.linspace(0,4,401)

	matplotlib.rc('xtick', labelsize=10)
	matplotlib.rc('ytick', labelsize=10)
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

	gs = gridspec.GridSpec(
        nrows=1, ncols=2, left=0.1, bottom=0.1, right=0.95, top=0.9,
        wspace=0.25, hspace=0., width_ratios=[1, 1])

	for i in range(len(rhofiles)):

                plt.figure(1)
                file = vxfiles[i]
                rhofile = rhofiles[i]

                n = file[8:]
                f = open(file, 'rb')
                rhof = open(rhofile, 'rb')

                X = f.read(num*num*size)
                X = np.array(struct.unpack(type*num*num, X))
                X = X.reshape((args.n, num))

                Y = rhof.read(num*num*size)
                Y = np.array(struct.unpack(type*num*num, Y))
                Y = Y.reshape((num,num))

                f.close()
                rhof.close()

                #plt.imshow(X)
                #plt.savefig('Plots/VXmap'+n+'.png')
                #print('Plots/VXmap'+n+'.png')

                #plt.cla()
                #plt.clf()

                #expected = np.cos(2*np.pi*np.linspace(0,1,1000).reshape(1000,1))/2*np.exp(-4*np.pi**2*visc*times[i])

                VX = X[:,0]/Y[:,0]
                y = np.linspace(1/2/len(VX), 1 - 1/2/len(VX), len(VX))
                ax1= plt.subplot(gs[0])
                #plt.plot(y, VX, 'k.--',label='1e-4')
                plt.imshow(X/Y,origin='lower',extent=(0,1,0,1))#,aspect='auto')
                #plt.legend(title='$\mu_r$')
                #plt.plot(np.linspace(0,1,1000), expected[:,0],'k')
                #plt.plot(np.linspace(0,1,1000), expected[:,1],'g')
                #plt.plot(np.linspace(0,1,1000), expected[:,2],'r')
                plt.xlabel('Position')
                #plt.ylabel('Momentum Density')
                plt.ylabel('Position')
                #plt.ylabel('Velocity')
                #plt.ylim(-0.5, 0.5)
                plt.title(f'Velocity (t = {times[i]:.2f}), '+r'$\mu_r =$ '+f'1e-4')
                plt.xticks([0,0.25,0.5,0.75,1])
                plt.yticks([0,0.25,0.5,0.75,1])
                divider = make_axes_locatable(plt.gca())
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(cax=cax,ticks=[-1/2,0,1/2])
                plt.clim(-0.5, 0.5)
                #plt.grid(alpha=0.2)

                #pdb.set_trace()
                plt.gcf().add_subplot(gs[1])#, sharey=ax1)
                #plt.subplot(gs[1])
                #plt.plot(y,Y[:,0],'k,--')
                plt.imshow(Y,origin='lower',extent=(0,1,0,1))#,aspect='auto')
                plt.xlabel('Position')
                plt.xticks([0,0.25,0.5,0.75,1])
                plt.yticks([0,0.25,0.5,0.75,1], ['','','','',''])
                plt.title(f'Density (t = {times[i]:.2f})')
                divider = make_axes_locatable(plt.gca())
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(cax=cax,ticks=[1,2,3])
                plt.clim(1,3)
                #plt.ylim(0,3)
                #plt.grid(alpha=0.2)

                plt.tight_layout()
                plt.savefig('Plots/VXprof'+n+'.png', dpi=230, bbox_inches = 'tight', pad_inches = 1/10)
                print('Plots/VXprof'+n+'.png')

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
		plt.savefig("Plots/phase"+n+".png")
		plt.figure()
		avg1 = F.mean(axis=0)
		avgnum = F.mean(axis=1)

		v = np.linspace(-10,10,num)
		plt.plot(v, avg2)
		plt.savefig('Plots/phaseavg'+n+'.png')
		print(f"Problem 6: Finished {n[1:]}")
		'''




if testproblem == 10:

	viscs = [-5, -3, -2, -1, 1]

	times = [0, 20, 41, 100]
	fig, axes = plt.subplots(len(times), 1, sharex=True, sharey=True, figsize=(4, 2*len(times)))
	for (ti, t) in enumerate(times):

		plt.figure("density", figsize=(12,8))
		plt.clf()

		for (vi, visc) in enumerate(viscs):

			file = f"Data_1e{visc}/rho_{t:04d}"
			f = open(file, 'rb')

			X = f.read(num*size)
			X = np.array(struct.unpack(type*num, X))

			axes[ti].set_xlim(0,1)
			axes[ti].set_ylim(0,15)
			axes[ti].plot(np.linspace(0.5/len(X),1-0.5/len(X), len(X)), X, label=rf"10$^{ {visc} }$")
			axes[ti].set_ylabel(f'Density (t ={t/200:.2f})', fontsize=12)
		axes[ti].grid()
		
	axes[0].legend(title=r"$\mu_r$")
	axes[-1].set_xlabel('x', fontsize=12)
	axes[0].set_title("Sine Wave")
	fig.subplots_adjust(hspace=0,wspace=1)
	fig.tight_layout()

	fig.savefig(f'Plots/Problem10/Rho.png', dpi=230)
	print(f'Problem 10: Plots/Rho.png')
	plt.close("fig")
	f.close()

	times = [0, 41, 120, 200]
	fig, axes = plt.subplots(len(viscs), len(times), sharex=True, sharey=True, figsize=(4*len(times), 4*len(viscs)))
	zmax = 0
	ims = []
	for (ti, t) in enumerate(times):

		axes[0][ti].set_title(f"time = {t/200:.2f}")
		for (vi, visc) in enumerate(viscs): 
		
			file = f"Data_1e{visc}/phase_{t:04d}"
			
			f = open(file, 'rb')

			F = f.read(num*numv*size)
			F = np.array(struct.unpack(type*num*numv,F))

			x = np.linspace(0, 1, num)
			v = np.linspace(-6, 6, numv)

			F = F.reshape((numv,num))
			Fmax = F.max()
			if Fmax > zmax:
				zmax = Fmax

			ims.append(axes[vi][ti].imshow(F, extent=[0,1,-6,6], origin="lower", aspect="auto", cmap="coolwarm"))
			axes[vi][ti].set_ylim(-2, 2)
			axes[vi][0].set_ylabel(rf"Velocity ($\mu_r = 10^{ {visc} }$)", fontsize=15)
		axes[vi][ti].set_xlabel("Position", fontsize=15)

	fig.subplots_adjust(hspace=0, wspace=0)

	# Make all axes have the same zlim
	for im in ims:
		im.set_clim(0, zmax/2)

	# add a colorbar for this axis
	cax = fig.add_axes([0.95, 0.05, 0.025, 0.9])
	fig.colorbar(im, cax=cax, orientation="vertical", cmap="coolwarm")

	fig.tight_layout(rect=(0,0,0.95,1))
	fig.savefig(f"Plots/Problem10/phase.png", dpi=230)
	print(f'Problem 10: Plots/phase.png')
	plt.close(fig)
		
		# for thing in things:
			# nullfmt = matplotlib.ticker.NullFormatter()         # no labels

			# # definitions for the axes
			# left, width = 0.15, 0.6
			# bottom, height = 0.15, 0.6
			# bottom_h = left_h = left + width + 0.05

			# rect_Phase = [left, bottom, width, height]
			# rect_avgx = [left, bottom_h, width, 0.15]
			# rect_avgy = [left_h, bottom, 0.15, height]

			# # start with a rectangular Figure
			# plt.figure(1, figsize=(12, 12))
			# plt.clf()

			# print(rect_Phase)
			# axPhase = plt.axes(rect_Phase)
			# axavgx = plt.axes(rect_avgx)
			# axavgy = plt.axes(rect_avgy)

			# # no labels
			# axavgx.xaxis.set_major_formatter(nullfmt)
			# axavgy.yaxis.set_major_formatter(nullfmt)

			# # the Phase plot:
			# axPhase.imshow(F, extent=[0,1,-6,6],aspect='auto')#,aspect=1/20)
			# axPhase.set_ylabel('Velocity')
			# axPhase.set_xlabel('Position')
			# axPhase.set_ylim((-6, 6))
			# axPhase.set_xlim((0,1))

			# #pdb.set_trace()
			# axavgx.plot(x, F.mean(axis=0))
			# axavgy.plot(F.mean(axis=1), v)#, orientation='horizontal')

			# axavgx.set_xlim(axPhase.get_xlim())
			# #axavgx.set_ylim(0,0.1)

			# #axavgy.set_xlim(0,1.5)
			# axavgy.set_ylim(axPhase.get_ylim())

			# plt.savefig(f"Plots/Problem10/phase_{i:04d}.png")
			# print(f"Problem 10: Finished {i:04d}")


if testproblem == 11:

	with imo.get_writer('Problem11.mp4', fps=30) as writer:

		
		for q, file in enumerate(rhofiles):

			fig, axes = plt.subplots(1,3, figsize=(15,4))

			buf = io.BytesIO()

			print(f"Working on file {q}")

			n = file[9:13]

			rho_file  = open(file, 'rb')
			rho = rho_file.read(num*size)
			rho = np.array(struct.unpack(type*num, rho))

			rhovx_file = open(f'Data/rhovx_{n}','rb')
			rho_vx = rhovx_file.read(num*size)
			rho_vx = np.array(struct.unpack(type*num, rho_vx))
			
			rhoE_file = open(f'Data/rhoE_{n}','rb')
			rho_E = rhoE_file.read(num*size)
			rho_E = np.array(struct.unpack(type*num, rho_E))

			axes[0].plot(np.linspace(0.5/len(rho),1-0.5/len(rho), len(rho)), rho)
			axes[0].set_xlim(0,1)
			axes[0].set_ylim(0,10)
			axes[0].set_xlabel('x')
			axes[0].set_ylabel('Density')
			axes[0].grid(alpha=0.5)

			vx = rho_vx/rho
			axes[1].plot(np.linspace(0.5/len(vx),1-0.5/len(vx), len(vx)), vx)
			axes[1].set_xlim(0,1)
			axes[1].set_ylim(0,10)
			axes[1].set_xlabel('x')
			axes[1].set_ylabel('Velocity')
			axes[1].grid(alpha=0.5)

			R = 1/2
			g = 7/5
			T = (g - 1)/R*(rho_E/rho - 0.5*rho_vx*rho_vx/rho/rho)
			print(T)
			axes[2].plot(np.linspace(0.5/len(T),1-0.5/len(T), len(T)), T)
			axes[2].set_xlim(0,1)
			axes[2].set_ylim(0,10)
			axes[2].set_xlabel('x')
			axes[2].set_ylabel('Temperature')
			axes[2].grid(alpha=0.5)
						
			plt.suptitle("Isothermal Shock")
			
			plt.savefig(buf,dpi=230)
			buf.seek(0)
			writer.append_data(imo.imread(buf))

			plt.cla()
			plt.clf()
			rho_file.close()
			rhovx_file.close()
			rhoE_file.close()
	
		writer.close()


if testproblem == 13:
    
	for file in rhofiles[:]:
		n = file[9:13]
		f = open(file, 'rb')

		X = f.read(num*size)
		X = np.array(struct.unpack(type*num, X))
		import pdb; pdb.set_trace()

		plt.xlim(0,1)
		plt.ylim(0,1.2)
		plt.plot(np.linspace(0.5/len(X),1-0.5/len(X), len(X)), X)#, 'o-')
		plt.grid()
		plt.savefig('Plots/Rho'+n+'.png')
		print('Plots/Rho'+n+'.png')

		plt.cla()
		plt.clf()
		f.close()