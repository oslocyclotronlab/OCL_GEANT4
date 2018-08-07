#!/usr/bin/python

import glob, os

multFactor = 1.15 # multiplication factor for cmp_bg

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
		print('Error: Creating directory. ' + directory)

outfolder = "../mama_spectra_smooth/"

createFolder(outfolder)

for file in glob.glob("cmp*"):
    os.system("cp %s %s/%s" % (file, outfolder, file))
    command = "./Smooth_and_Multiply.sh %s %s %f" % (outfolder, file, multFactor)
    os.system(command)

    # correct for some mistakes in writing the mama file
    command ="""vim -E -s {} <<-EOF 
:%s/\\n\\n!CAL/\\r!CAL/g
:x
EOF
""".format(outfolder+file)
    # print command
    os.system(command)