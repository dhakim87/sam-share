import glob
import os
import subprocess
import sys

if len(sys.argv) < 4:
    print("Usage: " + sys.argv[0] + " <in folder> <out folder> <executable>")
    exit(-1)

print("In Folder: " + sys.argv[1])
print("Out Folder: " + sys.argv[2])
print("Executable Path: " + sys.argv[3])

job_id = 0
num_jobs = 1

if len(sys.argv) == 6:
    job_id = int(sys.argv[4])
    num_jobs = int(sys.argv[5])

sams = glob.glob(sys.argv[1] + "/*.sam")
xzsams = glob.glob(sys.argv[1] + "/*.sam.xz")

all_sams = sams + xzsams

i = -1
for f in sorted(all_sams):
    i += 1
    if i % num_jobs != job_id:
        continue

    if f.endswith(".sam"):
        subprocess.run("cat " + f + " | " + sys.argv[3] + " > " + sys.argv[2] + "/" + os.path.basename(f)[:-4] + ".outsam", shell=True)
    elif f.endswith(".sam.xz"):
        subprocess.run("xzcat " + f + " | " + sys.argv[3] + " > " + sys.argv[2] + "/" + os.path.basename(f)[:-7] + ".outsam", shell=True)
