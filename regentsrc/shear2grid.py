import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('agg')
import matplotlib.pyplot as plt
import struct
import glob
import argparse
import pdb
#from sod import solve
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from collisionless import Collisionless
import h5py

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=int, default=0)
parser.add_argument('-n', type=int, default=128)
args = parser.parse_args()

size = 8
num  = args.n

type = 'd' #d is double, f is float, i is integer

testproblem = args.p

CBC = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

def FindFiles(Dir):
	rhofiles = glob.glob(Dir+'/rho_*')
	vxfiles = glob.glob(Dir+'/rhovx_*')
	vyfiles = glob.glob(Dir+'/rhovy_*')
	vzfiles = glob.glob(Dir+'/rhovz_*')
	Efiles = glob.glob(Dir+'/rhoE_*')

	rhofiles.sort()
	vxfiles.sort()
	vyfiles.sort()
	vzfiles.sort()
	Efiles.sort()

	#phasefiles = glob.glob(Dir+'/phase_*')
	#phasefiles.sort()

	return rhofiles, vxfiles, vyfiles, vzfiles, Efiles#, phasefiles

def Loadrho(rhofiles, i, size, num, type, df = 1):
	rho = rhofiles[i]
	rho = open(rho, 'rb')
	rho = rho.read(size*num*df)
	rho = np.array(struct.unpack(type*num*df, rho))[::df]
	return rho

rhofiles023, vxfiles023, vyfiles023, vzfiles023, Efiles023 = FindFiles('sh_mu1e-0_Pr2_3/Data')
rhofiles223, vxfiles223, vyfiles223, vzfiles223, Efiles223 = FindFiles('sh_mu1e-2_Pr2_3/Data')
rhofiles323, vxfiles323, vyfiles323, vzfiles323, Efiles323 = FindFiles('sh_mu1e-3_Pr2_3/Data')
rhofiles423, vxfiles423, vyfiles423, vzfiles423, Efiles423 = FindFiles('sh_mu1e-4_Pr2_3/Data')
rhofiles623, vxfiles623, vyfiles623, vzfiles623, Efiles623 = FindFiles('sh_mu1e-6_Pr2_3/Data')

rhofiles31, vxfiles31, vyfiles31, vzfiles31, Efiles31 = FindFiles('sh_mu1e-3_Pr1/Data')
rhofiles332, vxfiles332, vyfiles332, vzfiles332, Efiles332 = FindFiles('sh_mu1e-3_Pr3_2/Data')


rhomin = 1.5
rhomax = 1.5
'''
rhos = rhofiles21 + rhofiles31 + rhofiles41 + rhofiles61 + rhofiles323 + rhofiles332 + rhofiles423 + rhofiles432
for file in rhos:
	f = open(file, 'rb')

	X = f.read(num*num*size)
	X = np.array(struct.unpack(type*num*num, X))

	rhomin = np.minimum(X.min(),rhomin)
	rhomax = np.maximum(X.max(),rhomax)

print(f'Max and Min density are {rhomax:.3f} and {rhomin:.3f}')
'''

'''
gamma = 7/5
dustFrac = 0.0
npts = 1024
t = 0.15
left_state = (1.0, 1.0, 0)
right_state = (0.1, 0.125, 0.)

# left_state and right_state set pressure, density and u (velocity)
# geometry sets left boundary on 0., right boundary on 1 and initial
# position of the shock xi on 0.5
# t is the time evolution for which positions and states in tube should be
# calculated
# gamma denotes specific heat
# note that gamma and npts are default parameters (1.4 and 500) in solve
# function
positions, regions, values = solve(left_state=left_state,
    right_state=right_state, geometry=(0., 1., 0.5), t=t,
    gamma=gamma, npts=npts, dustFrac=dustFrac)
'''

viscs = [1e0, 1e-2, 1e-3, 1e-4, 1e-6]
labels = [0, -2, -3, -4, -6]
labels = [r'$10^{0}$', r'$10^{-2}$', r'$10^{-3}$', r'$10^{-4}$', r'$10^{-6}$']

viscs.reverse()
labels.reverse()
plt.figure(figsize=(10,4))

#Time Series
ts = [[] for i in range(len(viscs))]
ets = [[] for i in range(len(viscs))]
nsts = [[] for i in range(len(viscs))]
for i in range(401):
	rho0 = Loadrho(rhofiles023, i, size, num, type)
	rho2 = Loadrho(rhofiles223, i, size, num, type)
	rho3 = Loadrho(rhofiles323, i, size, num, type)
	rho4 = Loadrho(rhofiles423, i, size, num, type,df=1)
	rho6 = Loadrho(rhofiles623, i, size, num, type)

	V0 = Loadrho(vxfiles023, i, size, num, type)/rho0
	V2 = Loadrho(vxfiles223, i, size, num, type)/rho2
	V3 = Loadrho(vxfiles323, i, size, num, type)/rho3
	V4 = Loadrho(vxfiles423, i, size, num, type,df=1)/rho4
	V6 = Loadrho(vxfiles623, i, size, num, type)/rho6

	Y0 = Loadrho(vyfiles023, i, size, num, type)/rho0
	Y2 = Loadrho(vyfiles223, i, size, num, type)/rho2
	Y3 = Loadrho(vyfiles323, i, size, num, type)/rho3
	Y4 = Loadrho(vyfiles423, i, size, num, type,df=1)/rho4
	Y6 = Loadrho(vyfiles623, i, size, num, type)/rho6

	E0 = Loadrho(Efiles023, i, size, num, type)/rho0
	E2 = Loadrho(Efiles223, i, size, num, type)/rho2
	E3 = Loadrho(Efiles323, i, size, num, type)/rho3
	E4 = Loadrho(Efiles423, i, size, num, type,df=1)/rho4
	E6 = Loadrho(Efiles623, i, size, num, type)/rho6

	N = len(rho0)

	x = np.linspace(1/2/N, 1-1/2/N, N)
	plt.suptitle(f'Uniform Density Shear Viscosity Test (t={4*i/400:0.2f}')
	#plt.subplot(221)
	#plt.title('Density')
	#plt.plot(values['x'], values['rho'], 'k', lw=1, alpha=1, label='Euler')
	#plt.plot(cx, rho, 'k--', lw=1, alpha=1, label='CL')
	#plt.plot(x, rho6, CBC[1], lw=1, alpha=0.8, label=r'$10^{-6}$')
	#plt.plot(x, rho4, CBC[0], lw=1, alpha=0.8, label=r'$10^{-4}$')
	#plt.plot(x, rho3, CBC[2], lw=1, alpha=0.8, label=r'$10^{-3}$')
	#plt.plot(x, rho2, CBC[3], lw=1, alpha=0.8, label=r'$10^{-2}$')
	#plt.plot(x, rho0, CBC[4], lw=1, alpha=0.8, label=r'$10^0$')
	#plt.legend(title=r'$\mu_r$')
	#plt.
	#plt.grid(alpha=0.3)

	plt.subplot(121)
	plt.title('Velocity')
	R = 0.5
	g = 7/5
	#plt.plot(values['x'], values['u'], 'k', lw=1, alpha=1, label='Euler')
	#plt.plot(cx, u, 'k--', lw=1, alpha=1, label='CL')
	t = i*4./400
	print(f'time[{i}] = {t:0.5f}')
	ts[0].append(np.max(V0))
	ts[1].append(np.max(V2))
	ts[2].append(np.max(V3))
	ts[3].append(np.max(V4))
	ts[4].append(np.max(V6))

	for q in range(len(viscs)):
		expected = 0.5*np.cos(2*np.pi*x)*np.exp(-4*np.pi**2*viscs[q]*t)
		ets[q].append(np.max(expected))
		if q == 1:
			j = 0
		elif q == 0:
			j = 1
		else:
			j = q
		if viscs[q] in [1e0,1e-2,1e-6]:
			plt.plot(x, expected, CBC[j], linestyle='dashed',lw=1, alpha=0.8)#, label=labels[i])
	plt.plot(x, V6, CBC[1], lw=1, alpha=0.8, label=r'$10^{-6}$')
	#plt.plot(x, V4, CBC[0], lw=1, alpha=0.8, label=r'$10^{-4}$')
	#plt.plot(x, V3, CBC[2], lw=1, alpha=0.8, label=r'$10^{-3}$')
	plt.plot(x, V2, CBC[3], lw=1, alpha=0.8, label=r'$10^{-2}$')
	plt.plot(x, V0, CBC[4], lw=1, alpha=0.8, label=r'$10^0$')
	plt.legend(title=r'$\mu_r$')
	plt.grid(alpha=0.3)

	#plt.subplot(224)
	#plt.title('Temperature')
	T0 = (g-1)/R*(E0-V0**2/2-Y0**2/2)
	T2 = (g-1)/R*(E2-V2**2/2-Y2**2/2)
	T3 = (g-1)/R*(E3-V3**2/2-Y3**2/2)
	T4 = (g-1)/R*(E4-V4**2/2-Y4**2/2)
	T6 = (g-1)/R*(E6-V6**2/2-Y6**2/2)

	#plt.plot(x, T6, CBC[1], lw=1, alpha=0.8, label=r'$10^{-6}$')
	#plt.plot(x, T4, CBC[0], lw=1, alpha=0.8, label=r'$10^{-4}$')
	#plt.plot(x, T3, CBC[2], lw=1, alpha=0.8, label=r'$10^{-3}$')
	#plt.plot(x, T2, CBC[3], lw=1, alpha=0.8, label=r'$10^{-2}$')
	#plt.plot(x, T0, CBC[4], lw=1, alpha=0.8, label=r'$10^0$')
	#plt.legend(title=r'$\mu_r$')
	#plt.grid(alpha=0.3)
	#plt.xlabel('Position')

	plt.subplot(122)
	plt.title('Pressure')
	P0 = rho0*R*T0
	P2 = rho2*R*T2
	P3 = rho3*R*T3
	P4 = rho4*R*T4
	P6 = rho6*R*T6

	plt.plot(x, P6, CBC[1], lw=1, alpha=0.8, label=r'$10^{-6}$')
	plt.plot(x, P4, CBC[0], lw=1, alpha=0.8, label=r'$10^{-4}$')
	plt.plot(x, P3, CBC[2], lw=1, alpha=0.8, label=r'$10^{-3}$')
	plt.plot(x, P2, CBC[3], lw=1, alpha=0.8, label=r'$10^{-2}$')
	plt.plot(x, P0, CBC[4], lw=1, alpha=0.8, label=r'$10^0$')
	plt.legend(title=r'$\mu_r$')
	plt.grid(alpha=0.3)
	plt.xlabel('Position')

	plt.tight_layout(rect=[0,0,1,0.975])
	plt.savefig(f'Plots/rho{i:04d}.png',format="png",dpi=330)
	print(f'rho{i:04d}')

	plt.cla()
	plt.clf()

for i in range(401):
	file0 = h5py.File(f'../../../athena-public-version/sh128/Data0/from_array.cons.{i:05d}.athdf','r')
	file2 = h5py.File(f'../../../athena-public-version/sh128/Data2/from_array.cons.{i:05d}.athdf','r')
	file3 = h5py.File(f'../../../athena-public-version/sh128/Data3/from_array.cons.{i:05d}.athdf','r')
	file4 = h5py.File(f'../../../athena-public-version/sh128/Data4/from_array.cons.{i:05d}.athdf','r')
	file6 = h5py.File(f'../../../athena-public-version/sh128/Data6/from_array.cons.{i:05d}.athdf','r')

	nsrho0,nsrhoE0,nsrhovx0,nsrhovy0,nsrhovz0 = file0['cons']
	nsrho2,nsrhoE2,nsrhovx2,nsrhovy2,nsrhovz2 = file2['cons']
	nsrho3,nsrhoE3,nsrhovx3,nsrhovy3,nsrhovz3 = file3['cons']
	nsrho4,nsrhoE4,nsrhovx4,nsrhovy4,nsrhovz4 = file4['cons']
	nsrho6,nsrhoE6,nsrhovx6,nsrhovy6,nsrhovz6 = file6['cons']

	v0 = nsrhovy0/nsrho0
	v2 = nsrhovy2/nsrho2
	v3 = nsrhovy3/nsrho3
	v4 = nsrhovy4/nsrho4
	v6 = nsrhovy6/nsrho6

	nsts[0].append(v6.max())
	nsts[1].append(v4.max())
	nsts[2].append(v3.max())
	nsts[3].append(v2.max())
	nsts[4].append(v0.max())

#pdb.set_trace()
plt.figure()
#viscs/nus were reversed
ts.reverse()
for q in range(len(viscs)):
	plt.title(r'Residual of max $v_{x}$ w/ expected')
	if q == 1:
		j = 0
	elif q == 0:
		j = 1
	else:
		j = q
	#plt.semilogy(np.linspace(0,4,401), 0.5-np.array(ets[q]), CBC[j], linestyle='dashed')
	plt.plot(np.linspace(0,4,401), np.array(ets[q])-np.array(nsts[q]), CBC[j], linestyle='dashed')
	plt.plot(np.linspace(0,4,401), np.array(ets[q])-np.array(ts[q]), CBC[j], label = labels[q])

	plt.legend()
	plt.xlabel('Time')
	plt.ylabel('Maximum Velocity')
	plt.grid(alpha=0.5)
	plt.savefig(f'expresidual{q}.png')

	plt.clf()
	plt.cla()

plt.title(r'Evolution of max $v_{y}$')
for q in range(len(viscs)):
	if q == 1:
		j = 0
	elif q == 0:
		j = 1
	else:
		j = q
	#plt.semilogy(np.linspace(0,4,401), 0.5-np.array(ets[q]), CBC[j], linestyle='dashed')
	plt.plot(np.linspace(0,4,401), np.array(nsts[q]), CBC[j], linestyle='dashed')
	plt.plot(np.linspace(0,4,401), np.array(ts[q]), CBC[j], label = labels[q])
plt.legend()
plt.xlabel('Time')
plt.ylabel('Maximum Velocity')
plt.grid(alpha=0.5)
plt.savefig('TE.png')

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
                #plt.savefig('Check/VXmap'+n+'.png')
                #print('Check/VXmap'+n+'.png')

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
                plt.savefig('Check/VXprof'+n+'.png', dpi=230, bbox_inches = 'tight', pad_inches = 1/10)
                print('Check/VXprof'+n+'.png')

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
		avgnum = F.mean(axis=1)

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
