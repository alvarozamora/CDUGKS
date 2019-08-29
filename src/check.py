import numpy as np
import matplotlib.pyplot as plt
import glob
import struct


size = 8
num  = 128
type = 'd' #d is double, f is float,


problem = int(np.genfromtxt('Data/index.txt'))

print("Problem = ", problem)
x = np.genfromtxt('Data/x.txt')

plt.figure(1)
plt.figure(2,figsize = (10,8))
if problem == 1:
	for i in range(len(glob.glob("Data/rho*"))):
	
		num = str(i)
		while len(num) < 4:
			num = '0' + num

		file = 'Data/rho'+num+'.txt'
		rho = np.genfromtxt(file)
	
		plt.plot(x, rho)
		plt.ylim(0,1.2)
		plt.xlim(0,1)
		plt.grid()
		plt.title("Sod Shock Tube")
		plt.xlabel("Distance")
		plt.ylabel("Density")
		plt.savefig("Check/check"+num+".png")

		plt.cla()

		if i%10 == 0:
			print("Saved figure", i)

if problem == 2:
	for i in range(len(glob.glob("Data/rho*"))):
	
		num = str(i)
		while len(num) < 4:
			num = '0' + num

		file = 'Data/rho'+num+'.txt'
		rho = np.genfromtxt(file)
	
		n = np.round(np.sqrt(len(rho)))
		n = int(n)
		rho = rho.reshape((n,n))

		plt.figure(1)
		plt.imshow(rho)
		plt.title("KHI Density")
		plt.xlabel("Distance (x)")
		plt.ylabel("Distance (y)")
		plt.savefig("Check/check"+num+".png")

		plt.cla()
		
		plt.figure(2)
		plt.plot(rho[:,4])
		plt.ylim(0.5, 2.5)
		plt.grid()
		plt.title("KHI Line")
		plt.xlabel("Distance")
		plt.ylabel("Density")
		plt.savefig("Check2/check"+num+".png")

		if i%10 == 0:
			print("Saved figure", i)

		plt.cla()

