import numpy as np
import h5py
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pymorton as pm
import struct

def loaddata(i):
	rho = open(f'Data/rho_{i:04d}','rb')
	vx = open(f'Data/rhovx_{i:04d}','rb')
	vy = open(f'Data/rhovy_{i:04d}','rb')
	E = open(f'Data/rhoE_{i:04d}','rb')

	rho = rho.read(4*2*100*100)
	vx = vx.read(4*2*100*100)
	vy = vy.read(4*2*100*100)
	E = E.read(4*2*100*100)

	rho = np.array(struct.unpack('f'*2*100*100,rho))[::2].reshape(100,100)
	vx = np.array(struct.unpack('f'*2*100*100,vx))[::2].reshape(100,100)/rho
	vy = np.array(struct.unpack('f'*2*100*100,vy))[::2].reshape(100,100)/rho
	E = np.array(struct.unpack('f'*2*100*100,E))[::2].reshape(100,100)

	return rho, vx, vy, E

plt.figure(figsize=(10,10))
dt = 2*np.pi/5/200
for i in range(200):

	rho, vx, vy, E = loaddata(i)
	x = y = np.linspace(1/100/2, 1-1/100/2, 100)-1/2
	xx, yy = np.meshgrid(x, y)
	rr = np.sqrt(xx**2 + yy**2)
	pp = np.arctan2(yy,xx)

	vr =  vx*np.cos(pp) + vy*np.sin(pp)
	vp = -vx*np.sin(pp) + vy*np.cos(pp)

	plt.subplot(221)
	plt.imshow(rho, origin='lower',extent = (0,1,0,1))
	plt.colorbar()
	plt.clim(0.9,1.1)
	plt.title('Density')

	plt.subplot(222)
	plt.imshow(vr, origin='lower',extent = (0,1,0,1))
	plt.colorbar()
	plt.clim(-1,1)
	plt.title('R-Velocity')

	plt.subplot(223)
	plt.imshow(vp,origin='lower',extent = (0,1,0,1))
	plt.colorbar()
	plt.clim(-1,1)
	plt.xlabel('Position')
	plt.title('Phi-Velocity')

	plt.subplot(224)
	plt.imshow(E,origin='lower',extent = (0,1,0,1))
	plt.colorbar()
	plt.title('Energy Density')

	plt.tight_layout(rect=[0,0,1,0.95])
	plt.suptitle(f"Gresho Vortex Problem (t = {dt*i:03.3f}, angle = {i*2*np.pi/200:.4f})")
	plt.savefig(f'Plots/2dpanel{i:05d}.png')
	print(f"Saved panel{i:05d}")
	plt.clf()
	plt.cla()
