import imageio as io
import os
import glob
import re
import numpy as np

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

print("Making Gif")
file_names = glob.glob("Plots/Rho*.png")
phase_names = glob.glob("Plots/phase_*.png")
sort_nicely(file_names)
sort_nicely(phase_names)
#print(file_names)
#making animation

problem = int(input("which problem (int)"))
phase_bool = False
if problem in [10]:
	phase_bool = bool(input("phase space (bool)?"))
#with io.get_writer('cdugks.gif', duration=dur) as writer:
#    for filename in file_names:
#        image = io.imread(filename)
#        writer.append_data(image)
#    writer.close()
filename = f"cdugks{problem}.mp4"
pfilename = f"cdugks{problem}_phase.mp4"
print(filename)
with io.get_writer(filename, fps=30) as writer:
	for fname in file_names:
		image = io.imread(fname)
		writer.append_data(image)
	writer.close()
print("MP4 Complete")

if phase_bool:
	with io.get_writer(pfilename, fps=30) as writer:
		for fname in phase_names:
			image = io.imread(fname)
			writer.append_data(image)
		writer.close()
	print("MP4 Complete")


