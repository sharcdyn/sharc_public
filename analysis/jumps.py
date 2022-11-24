#!/bin/env python3
import numpy
import os

jumps_from=int(input("From which potential are the jumps (0) all "))
jumps_to=int(input("To which potential are the jumps (0) all "))

fjump=open("jumps.xyz","w")
for i in range(1000):
 fich="TRAJ"+str(i)+"/traj.xyz"
 if os.path.isfile(fich):
  pots=[]
  comments=[]
  geoms=[]
  with open(fich) as f:
   while True:
    line=f.readline()
    kk=line.split()
    if len(kk) == 0:
     break
    nat=int(kk[0])
    line=f.readline()
    kk=line.split()
    pots.append(int(kk[3]))
    comments.append(line)
    xyz=[]
    for j in range(nat):
     xyz.append(f.readline())
    geoms.append(xyz)
  print(pots)
  for j in range(1,len(pots)):
   if pots[j-1] != pots[j]:
    if pots[j-1]==jumps_from or jumps_from==0:
     if pots[j]==jumps_to or jumps_to==0:
      fjump.write(str(nat)+"\n")
      fjump.write("TRAJ"+str(i)+" from "+str(pots[j-1])+" "+comments[j])
      line=""
      for k in geoms[j]:
       line+=k
      fjump.write(line)
    
  
