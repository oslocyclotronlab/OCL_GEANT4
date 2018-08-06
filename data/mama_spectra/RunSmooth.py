#!/usr/bin/python

import glob, os

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
print ('Error: Creating directory. ' + directory)

outfolder = "../mama_spectra_smooth/"

createFolder(outfolder)

for file in glob.glob("cmp*"):
    os.system("cp %s %s/%s" % (file, outfolder, file))
    command = "./Smooth.sh %s %s" % (outfolder, file)
    os.system(command)