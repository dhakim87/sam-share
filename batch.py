import glob
import os
import subprocess
import sys

if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <in folder> <out folder> <executable>")

print("In Folder: " + sys.argv[1])
print("Out Folder: " + sys.argv[2])
print("Executable Path: " + sys.argv[3])

sams = glob.glob(sys.argv[1] + "/*.sam")
xzsams = glob.glob(sys.argv[1] + "/*.sam.xz")

all_sams = sams + xzsams

for f in all_sams:
    if f.endswith(".sam"):
        subprocess.run("echo cat " + f + " | " + sys.argv[3] + " > " + sys.argv[2] + "/" + os.path.basename(f)[:-4] + ".outsam", shell=True)
    elif f.endswith(".sam.xz"):
        subprocess.run("echo xzcat " + f + " | " + sys.argv[3] + " > " + sys.argv[2] + "/" + os.path.basename(f)[:-7] + ".outsam", shell=True)
