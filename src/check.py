import numpy as np
import matplotlib.pyplot as plt
import glob
import struct


size = 8
num  = 128
type = 'd' #d is double, f is float,


x = np.genfromtxt('Data/x.txt')

plt.figure()
plt.grid()
plt.xlabel("Distance")
plt.ylabel("Density")
plt.title("Sod Shock Tube")


for i in range(len(glob.glob("Data/rho*"))):
	
	num = str(i)
	while len(num) < 4:
		num = '0' + num

	file = 'Data/rho'+num+'.txt'
	rho = np.genfromtxt(file)
	
	plt.plot(x, rho)
	plt.savefig("Check/check"+num+".png")

	plt.cla()

	if i%10 == 0:
		print("Saved figure", i)
