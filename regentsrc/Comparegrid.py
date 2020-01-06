import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import struct
import glob
import pdb

files = glob.glob('Data/rho*')
files2 = glob.glob('Data2/rho*')
size = 8
num  = 256
type = 'd' #d is double, f is float, i is integer

testproblem = 1

files.sort()
files2.sort()
plt.figure()

file = files[-1]
file2 = files2[-1]
f = open(file, 'rb')
f2 = open(file2, 'rb')

X = f.read(num*size)
X = np.array(struct.unpack(type*num, X))
X2 = f2.read(num*size)
X2 = np.array(struct.unpack(type*num, X2))
plt.title("Residual")
plt.semilogy(X-X2)
plt.savefig("Residual.png")

pdb.set_trace()
if testproblem == 1:
	for i in range(len(files)):

                file = files[i]
                file2 = files2[i]

                n = file[10:14]
                f = open(file, 'rb')
                f2 = open(file2, 'rb')

                X = f.read(num*size)
                X = np.array(struct.unpack(type*num, X))
                X2 = f2.read(num*size)
                X2 = np.array(struct.unpack(type*num, X2))

                plt.title("Time = "+str(np.round((i+1)/len(files),3)))
                plt.xlim(0,1)
                plt.ylim(0,1.2)
                plt.plot(np.linspace(0.5/len(X), 1-0.5/len(X),len(X)),X)
                plt.plot(np.linspace(0.5/len(X2), 1-0.5/len(X2),len(X2)), X2)
                plt.grid()
                plt.savefig('Galilean/Rho'+n+'.png')
                print('Galilean/Rho'+n+'.png')

                plt.cla()
                plt.clf()
                f.close()
                f2.close()

		
