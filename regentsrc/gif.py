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
file_names = glob.glob("Check/Rho*.png")
phase_names = glob.glob("Check/phase_*.png")
sort_nicely(file_names)
sort_nicely(phase_names)
#print(file_names)
#making animation

import pdb
pdb.set_trace()

if problem != 2:
	dur = 1/30.
elif problem == 2:
	dur = 1/30.
	file_names2 = glob.glob("Check2/*.png")
	sort_nicely(file_names2)

#with io.get_writer('cdugks.gif', duration=dur) as writer:
#    for filename in file_names:
#        image = io.imread(filename)
#        writer.append_data(image)
#    writer.close()


if problem == 2:
	with io.get_writer('cdugks2.gif', duration=dur) as writer:
		for filename in file_names2:
			image = io.imread(filename)
			writer.append_data(image)
		writer.close()
if problem == 6:
	with io.get_writer('cdugks6.gif', duration=dur) as writer:
		for filename in phase_names:
			image = io.imread(filename)
			writer.append_data(image)
		writer.close()
	print("GIF Complete")
	with io.get_writer('cdugks6.mp4', duration=dur) as writer:
		for filename in phase_names:
			image = io.imread(filename)
			writer.append_data(image)
		writer.close()
	print("MP4 Complete")

