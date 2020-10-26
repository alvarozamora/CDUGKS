import numpy as np
import h5py
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
rhofiles123, vxfiles123, vyfiles123, vzfiles123, Efiles123 = FindFiles('sh_mu1e-1_Pr2_3/Data')
rhofiles223, vxfiles223, vyfiles223, vzfiles223, Efiles223 = FindFiles('sh_mu1e-2_Pr2_3/Data')
rhofiles323, vxfiles323, vyfiles323, vzfiles323, Efiles323 = FindFiles('sh_mu1e-3_Pr2_3/Data')

viscs = [1e0, 1e-1, 1e-2,1e-3]
labels = [r'$10^{0}$', r'$10^{-1}$', r'$10^{-2}$', r'$10^{-3}$']

plt.figure()
#Time Series
ts = [[] for i in range(len(viscs))]
ets = [[] for i in range(len(viscs))]
nsts = [[] for i in range(len(viscs))]
for i in range(401):
	file0 = h5py.File(f'../../../athena-public-version/sh128/Data0/from_array.cons.{i:05d}.athdf','r')
	file1 = h5py.File(f'../../../athena-public-version/sh128/Data1/from_array.cons.{i:05d}.athdf','r')
	file2 = h5py.File(f'../../../athena-public-version/sh128/Data2/from_array.cons.{i:05d}.athdf','r')
	file3 = h5py.File(f'../../../athena-public-version/sh128/Data3/from_array.cons.{i:05d}.athdf','r')

	nsrho0,nsrhoE0,nsrhovx0,nsrhovy0,nsrhovz0 = file0['cons']
	nsrho1,nsrhoE1,nsrhovx1,nsrhovy1,nsrhovz1 = file1['cons']
	nsrho2,nsrhoE2,nsrhovx2,nsrhovy2,nsrhovz2 = file2['cons']
	nsrho3,nsrhoE3,nsrhovx3,nsrhovy3,nsrhovz3 = file3['cons']

	v0 = nsrhovy0/nsrho0
	v1 = nsrhovy1/nsrho1
	v2 = nsrhovy2/nsrho2
	v3 = nsrhovy3/nsrho3

	nsts[0].append(v0.max())
	nsts[1].append(v1.max())
	nsts[2].append(v2.max())
	nsts[3].append(v3.max())

for i in range(401):
	rho0 = Loadrho(rhofiles023, i, size, num, type)
	rho1 = Loadrho(rhofiles123, i, size, num, type)
	rho2 = Loadrho(rhofiles223, i, size, num, type)
	rho3 = Loadrho(rhofiles323, i, size, num, type)

	V0 = Loadrho(vxfiles023, i, size, num, type)/rho0
	V1 = Loadrho(vxfiles123, i, size, num, type)/rho1
	V2 = Loadrho(vxfiles223, i, size, num, type)/rho2
	V3 = Loadrho(vxfiles323, i, size, num, type)/rho3

	Y0 = Loadrho(vyfiles023, i, size, num, type)/rho0
	Y1 = Loadrho(vyfiles123, i, size, num, type)/rho1
	Y2 = Loadrho(vyfiles223, i, size, num, type)/rho2
	Y3 = Loadrho(vyfiles323, i, size, num, type)/rho3

	E0 = Loadrho(Efiles023, i, size, num, type)
	E1 = Loadrho(Efiles123, i, size, num, type)
	E2 = Loadrho(Efiles223, i, size, num, type)
	E3 = Loadrho(Efiles323, i, size, num, type)
	N = len(rho0)

	ts[0].append(np.max(V0))
	ts[1].append(np.max(V1))
	ts[2].append(np.max(V2))
	ts[3].append(np.max(V3))

	t = i*4.0/400
	for q in range(len(viscs)):
		expected = 0.5*np.exp(-4*np.pi**2*viscs[q]*t)
		ets[q].append(expected)
plt.title(r'Evolution of max $v_{y}$')
pdb.set_trace()
for q in range(len(viscs)):
	plt.plot(np.linspace(0,4,401), np.array(nsts[q]), CBC[q], linestyle='dashed')
	plt.plot(np.linspace(0,4,401), np.array(ts[q]), CBC[q], label = labels[q])
plt.legend(title=r'$\mu_r$',loc=0)
plt.xlabel('Time')
plt.xlim(0,2)
plt.ylim(0,0.51)
plt.ylabel('Maximum Velocity')
plt.grid(alpha=0.5)
plt.savefig('TE.png')
