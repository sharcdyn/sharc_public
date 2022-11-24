import sys
import numpy as np
import os

listoftrajs=[]
if len(sys.argv)>1:
 kk=np.loadtxt(sys.argv[1])
 for n in kk:
  listoftrajs.append(int(n))
else:
 n=0
 while True:
  fich="../TRAJ"+str(n)+"/traj.xyz"
  if os.path.isfile(fich):
   listoftrajs.append(n)
   n+=1
  else:
   if n!=0:
    break

print(listoftrajs)

Ekins=[]
Epots=[]
Etots=[]
for n in listoftrajs:
 fich="../TRAJ"+str(n)+"/traj.xyz"
 print("Doing ",fich)
 os.system("ln -sf "+fich+" .")
 os.system("./sharc_ovnm.exe >/dev/null")
 kk=np.loadtxt("Ekin_harmonic.out")
 print(kk.shape)
 if n==0:
  maxsteps=kk.shape[0]
  times=kk[:,0]
  nvib=kk.shape[1]-1
  Ek_av=np.zeros(nvib)
  nforstep=np.zeros( maxsteps )
  Ek_av=np.average(kk,axis=0)
  kk=np.loadtxt("Epot_harmonic.out")
  Ep_av=np.average(kk,axis=0)
  kk=np.loadtxt("Etot_harmonic.out")
  Et_av=np.average(kk,axis=0)
  np.savetxt("Ekin0_t.dat",Ek_av)
  np.savetxt("Epot0_t.dat",Ep_av)
  np.savetxt("Etot0_t.dat",Et_av)
  os.system("mv Ekin_harmonic.out Ekin0.dat")
  os.system("mv Epot_harmonic.out Epot0.dat")
  os.system("mv Etot_harmonic.out Etot0.dat")
 else:
  if kk.shape[0]>maxsteps:
   newmaxsteps=kk.shape[0]
   mforstep=np.zeros( newmaxsteps )
   mforstep[:maxsteps]=nforstep[:]
   times=kk[:,0]*1.
   nforstep=mforstep*1.
   maxsteps=newmaxsteps*1
    
  for i in range(kk.shape[0]):
   nforstep[i]+=1
  Ekins.append(kk)
  Epots.append(np.loadtxt("Epot_harmonic.out"))
  Etots.append(np.loadtxt("Etot_harmonic.out"))

print(len(Ekins),maxsteps,nforstep)

Ekin_av=np.zeros( (maxsteps,nvib+1) )
Ekin_av[:,0]=times
Epot_av=np.zeros( (maxsteps,nvib+1) )
Epot_av[:,0]=times
Etot_av=np.zeros( (maxsteps,nvib+1) )
Etot_av[:,0]=times
for i in range(maxsteps):
 for (j,k,l) in zip(Ekins,Epots,Etots):
  if j.shape[0]>i:
   for n in range(nvib):
    Ekin_av[i,n+1]+=j[i,n+1]/nforstep[i]
    Epot_av[i,n+1]+=k[i,n+1]/nforstep[i]
    Etot_av[i,n+1]+=l[i,n+1]/nforstep[i]
  
np.savetxt("Ekin_average.dat",Ekin_av)
Ek_av=np.average(Ekin_av,axis=0)
np.savetxt("Ekin_av_t.dat",Ek_av)

np.savetxt("Epot_average.dat",Epot_av)
Ep_av=np.average(Epot_av,axis=0)
np.savetxt("Epot_av_t.dat",Ep_av)

np.savetxt("Etot_average.dat",Etot_av)
Et_av=np.average(Etot_av,axis=0)
np.savetxt("Etot_av_t.dat",Et_av)
