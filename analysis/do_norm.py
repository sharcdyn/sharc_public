#!/usr/bin/env python3

import numpy as np
import os,sys

dE=1./27.211
Erest=False

traj0=int(sys.argv[1])
trajf=int(sys.argv[2])
if len(sys.argv)>3:
 dE=float(sys.argv[3])/27.211

if len(sys.argv)>4:
 Erest=True
 Emin=float(sys.argv[4])/27.211
 Emax=float(sys.argv[5])/27.211


print(sys.argv)
print(traj0,trajf,dE)
if Erest:
 print(Emin,Emax)

cont_traj=False
if trajf<0:
 cont_traj=True
 trajf=-trajf

trajs=[]
valids=[]
invalids=[]
ntrajs=0
for itraj in range(traj0,trajf+1):
 print(itraj)
 fich="TRAJ"+str(itraj)+"/traj.xyz"
 if os.path.isfile(fich):
  with open(fich) as f:
   traj={}
   traj["traj"]=itraj
   traj["geom"]=[]
   traj["vel"]=[]
   traj["title"]=[]
   traj["sat"]=[]
   traj["pot"]=[]
   traj["time"]=[]
   geom=[]
   vel=[]
   sat=[]
   n=-2
   for line in f:
    n+=1
    kk=line.replace("\n","").split()
    if n==-1:
     nat=int(kk[0])
    if n==0:
     title=line
     traj["title"].append(title)
     traj["pot"].append(int(kk[3]))
     traj["time"].append(float(kk[2]))
    if n>=1 and n<=nat:
     sat.append(kk[0])
     for i in range(3):
      geom.append(float(kk[i+1]))
      vel.append(float(kk[i+4])/1000.)
    if n==nat:
     n=-2
     traj["geom"].append(geom)
     traj["vel"].append(vel)
     traj["sat"].append(sat)
     geom=[]
     vel=[]
     sat=[]

  valid=True
  with open("TRAJ"+str(itraj)+"/coef.out") as f:
   for line in f:
    kk=line.split()
    if np.isnan(float(kk[1])):
     valid=False
  with open("TRAJ"+str(itraj)+"/energy.out") as f:
   energy0=None
   for line in f:
    kk=line.split()
    energy=float(kk[1])+float(kk[2])
    if energy0 is None:
     energy0=energy
     venergy=float(kk[2])-float(kk[3])
    if abs(energy-energy0)>dE:
     valid=False
    if Erest:
     if venergy<Emin or venergy>Emax:
      valid=False
 
  for files in ["pop_raw","pop_a"]:
   lasttime=len(traj["time"])
   with open("TRAJ"+str(itraj)+"/"+files+".out") as f:
    traj[files]=[]
    for line in f:
     pop=[]
     kk=line.replace("\n","").split()
     for i in range(1,len(kk)):
      if np.isnan(float(kk[i])):
       valid=False
      pop.append(float(kk[i]))
###
     if len(traj[files])<=lasttime-1:
      if abs(float(kk[0])-traj["time"][len(traj[files])])<=1e-6:
       traj[files].append(pop.copy())
  
  if valid:
   valids.append(ntrajs)
  else:
   invalids.append(itraj)
  trajs.append(traj.copy())
  ntrajs+=1

maxsteps=0
maxpot=0
imax=0
for i in valids:
 if maxsteps <= len(trajs[i]["sat"]):
  maxsteps=len(trajs[i]["sat"])
  imax=i
 for j in trajs[i]["pot"]:
  if maxpot <= j:
   maxpot=j
 for j in trajs[i]:
  if j != "traj":
   print(i,j,len(trajs[i][j]))
  
times=trajs[imax]["time"]

kk=[]
for i in valids:
 kk.append(trajs[i]["traj"])
print("Valid trajs in file valids.dat",kk)
np.savetxt("valids.dat",kk)
print("Total trajs ",ntrajs)
print("Invalid trajs ",invalids)
print("Number of valid trajs ",len(valids))
print(maxsteps,maxpot)

if not os.path.isdir("steps"):
 os.mkdir("steps")
else:
 filestodelete=os.listdir("steps")
 for i in filestodelete:
  if "xyz" in i:
   os.remove("steps/"+i)
if cont_traj:
 print("Continuing the trajectory using the last point")
 files=["pot","pop_raw","pop_a"]
 for i in files:
  for itraj in valids:
   maxstep=len(trajs[itraj][i])
   if i=="pot":
    kk=trajs[itraj][i][-1]
   else:
    kk=trajs[itraj][i][-1].copy()
   for istep in range(maxstep,maxsteps):
    trajs[itraj][i].append(kk)

f1=open("steps/norm.out","w")
f2=open("steps/norm_raw.out","w")
f3=open("steps/norm_a.out","w")
for istep in range(maxsteps):
 norm=np.zeros( maxpot )
 norm_raw=np.zeros( len(trajs[0]["pop_raw"][0]) )
 norm_a=np.zeros( len(trajs[0]["pop_a"][0]) )
 nnorm=0
 nnorm_a=0
 nnorm_raw=0

 for itraj in valids:
  if len(trajs[itraj]["pot"])>istep:
   nnorm+=1
   pot=trajs[itraj]["pot"][istep]-1
   norm[int(pot)]+=1.

  if len(trajs[itraj]["pop_raw"])>istep:
   nnorm_raw+=1
   for i in range(len(norm)):
    norm_raw[i]+=trajs[itraj]["pop_raw"][istep][i]*1.

  if len(trajs[itraj]["pop_a"])>istep:
   nnorm_a+=1
   for i in range(len(norm)):
    norm_a[i]+=trajs[itraj]["pop_a"][istep][i]*1.
  
# WRITING norm
 f1.write(str(times[istep])+" "+str(nnorm))
 for i in range(maxpot):
  f1.write(" "+str(norm[i]/float(nnorm)) )
 f1.write("\n")

# WRITING RAW norm
 f2.write(str(times[istep]))
 for i in norm_raw:
  f2.write(" "+str(i/float(nnorm_raw)) )
 f2.write("\n")

# WRITING ADIABATIC norm
 f3.write(str(times[istep]))
 for i in norm_a:
  f3.write(" "+str(i/float(nnorm_a)) )
 f3.write("\n")

 g={}
 for ipot in range(maxpot):
  if norm[ipot]>0:
   fich="steps/step{:05d}-{:03d}.xyz".format(istep,ipot+1)
   g[ipot]=open(fich,"w")
 for itraj in valids:
  if len(trajs[itraj]["title"])>istep:
   pot=trajs[itraj]["pot"][istep]-1
   for ipot in range(maxpot):
    if pot==ipot:
     f=g[pot]
     f.write(str(nat)+"\n")
     kk=trajs[itraj]["title"][istep].split()

# Changing the comment line to include the traj number instead time (2nd and 3rd fields)
     towrite=""
     for i in range(len(kk)):
      if i==1:
       towrite+="TRAJ "
      elif i==2:
       towrite+=str(trajs[itraj]["traj"])+" "
      elif i==3:
       towrite+=""
      elif i==len(kk)-1:
       towrite+=kk[i]+"\n"
      else :
       towrite+=kk[i]+" "

     f.write(towrite)
     sat=trajs[itraj]["sat"][istep]
     geom=trajs[itraj]["geom"][istep]
     vel=trajs[itraj]["vel"][istep]
     n=0
     for iat in range(len(sat)):
      kk=sat[iat]
      for i in range(3):
       kk+=" "
       kk+=str(geom[n+i])
      for i in range(3):
       kk+=" "
       kk+=str(vel[n+i]*1000.)
      n+=3
      kk+="\n"
      f.write(kk)
 for ipot in g:
  g[ipot].close() 

f1.close()
f2.close()
f3.close()
