import sys,os
import numpy as np

Ew=float(sys.argv[1])/27.211
valids=[]
if os.path.isfile("valids.dat"):
 kk=np.loadtxt("valids.dat")
 if (len(kk.shape)==0):
  valids.append(int(kk))
 else:
  for i in kk:
   valids.append(int(i))
else:
 traj0=int(sys.argv[2])
 trajf=int(sys.argv[3])
 for i in range(traj0,trajf+1):
  valids.append(i)

print(valids)

def read_dip(fich):
 kk=np.loadtxt(fich)
 ntimes0=kk.shape[0]
 npot0=int( (kk.shape[1]-1)/2 )
 dip=np.zeros( ( ntimes0,npot0 ), dtype=np.complex )
 for i in range(ntimes0):
  for j in range(npot0):
   dip[i,j]=kk[i,2*j+1]+1j*kk[i,2*j+2]
 npot1=int( np.sqrt(npot0) )
 dip=dip.reshape( ( ntimes,npot1,npot1) )
 return dip

def smooth(a,w,b):
 import scipy.signal
 dx=(b[len(b)-1]-b[0])/(len(b)-1)
 width=2*np.sqrt(2*np.log(2))
 width=w/width /dx
 kernel=scipy.signal.gaussian(len(a),width)
 data=np.convolve(a,kernel,"same") #/sum(kernel)
 return data

def assign(grid_x,grid_y,prob=None):
 new_grid=np.zeros(len(grid_x))
 if prob.all()==None:
  kk=np.ones(len(grid_y))
 else:
  kk=prob
       
 for x in range(len(grid_y)):
  idx=(np.abs(grid_x-grid_y[x]).argmin())
  new_grid[idx] +=kk[x]
 return new_grid

npot=None
trajs=[]
times=[]
for itraj in valids:
 traj={}
 dir0="TRAJ"+str(itraj)+"/"
 if os.path.isdir(dir0):
  kk=np.loadtxt(dir0+"energy.out")
  print("Doing TRAJ ",itraj,kk.shape)
  time=kk[:,0]
  ntimes=time.shape[0]
  energy=kk[:,3:]
  if ntimes>len(times):
   times=time*1.
  if npot==None:
   npot=energy.shape[1]
  else:
   if npot!=energy.shape[1]:
    print("The number of potentials is not the same in all trajs")
    sys.exit()

  Ek=kk[:,1]
  pot=[]
  for i in range(ntimes):
   for j in range(npot):
    if kk[i,2]==energy[i,j]:
     pot.append(j)

  dip1=read_dip(dir0+"dipole1.out")
  dip2=read_dip(dir0+"dipole2.out")
  dip3=read_dip(dir0+"dipole3.out")
  traj['traj']=int(itraj)
  traj['time']=time
  traj['energy']=energy
  traj['dip1']=dip1
  traj['dip2']=dip2
  traj['dip3']=dip3
  traj['Ek']=Ek
  traj['pot']=pot
  trajs.append(traj)

print("Data ready")
ntimes=len(times)
print("Times in the longest trajectory ",ntimes,times[0],times[1],"...",times[-2],times[-1])

oscs=[]
dEs=[]
dEs_tot=[]
for i in range(ntimes):
 oscs_time=[]
 dEs_time=[]
 ntraj=0
 for j in trajs:
  time=j['time']
  if time[-1] >= times[i]:
   ntraj+=1
   energy=j['energy'][i,:]
   pot=j['pot'][i]
   dip=np.zeros( (3,npot,npot), dtype=np.complex )
   dip[0,:,:]=j['dip1'][i,:,:]
   dip[1,:,:]=j['dip2'][i,:,:]
   dip[2,:,:]=j['dip3'][i,:,:]
   dE=np.zeros( npot )
   osc=np.zeros( (3,npot) )
   for k in range(npot):
    dE[k]=energy[k]-energy[pot]
    for l in range(3):
     osc[l,k]=np.abs(dip[l,k,pot])**2
     osc[l,k]=osc[l,k]*dE[k]*2./3.
   dEs_time.append(dE)
   oscs_time.append(osc)
   for iE in dE:
    dEs_tot.append(iE)
 dEs.append(dEs_time)
 oscs.append(oscs_time)

del trajs

print("Creating energy grid")
Emin=np.min(np.array(dEs_tot))
print("Minimum energy ",Emin*27.211,(Emin-10*Ew)*27.211)
Emax=np.max(np.array(dEs_tot))
print("Maximum energy ",Emax*27.211,(Emax+10*Ew)*27.211)

gridE=np.arange(Emin-10*Ew,Emax+10*Ew,Ew/10.)
print("Energy grid ",np.min(gridE),np.max(gridE),len(gridE))

probs=np.zeros( (ntimes,3,npot,len(gridE)) )
for itime in range(ntimes):
 dE=np.array(dEs[itime])
 osc=np.array(oscs[itime])
 for k in range(npot):
  for l in range(3):
   gridEs=assign(gridE,dE[:,k], np.abs(osc[:,l,k]) )
   probs[itime,l,k,:]=smooth(gridEs,Ew,gridE)

import gzip
for l in range(4):
 filename="td-spectra"+str(l+1)+".out"
 if l==3:
  filename="td-spectra.out"
 print("Saving ",filename+".gz")
 with gzip.open(filename+".gz","wt",compresslevel=9) as f:
  spec=np.zeros( (npot+3,len(gridE)) ) 
  spec[1,:]=gridE
  for itime in range(ntimes):
   spec[0,:]=times[itime]*1.
   if l==3:
    spec[3:,:]=np.sum(probs[itime,:,:,:],axis=0)
   else:
    spec[3:,:]=probs[itime,l,:,:]
   spec[2,:]=np.sum(spec[3:,:],axis=0)
   np.savetxt(f,np.transpose(spec))
   f.write("\n")


print("Creating experimental energy grid")
Emin=0.
Emax=np.max(np.abs(np.array(dEs_tot)))
print("Maximum energy ",Emax*27.211,(Emax+10*Ew)*27.211)

gridE=np.arange(Emin,Emax+10*Ew,Ew/10.)
print("Energy grid ",np.min(gridE),np.max(gridE),len(gridE))

probs=np.zeros( (ntimes,3,npot,len(gridE)) )
for itime in range(ntimes):
 dE=np.array(dEs[itime])
 osc=np.array(oscs[itime])
 for k in range(npot):
  for l in range(3):
   gridEs=assign(gridE,np.abs(dE[:,k]), np.abs(osc[:,l,k]) )
   probs[itime,l,k,:]=smooth(gridEs,Ew,gridE)

import gzip
for l in range(4):
 filename="td-spectra"+str(l+1)+"_exp.out"
 if l==3:
  filename="td-spectra_exp.out"
 print("Saving ",filename+".gz")
 with gzip.open(filename+".gz","wt",compresslevel=9) as f:
  spec=np.zeros( (npot+3,len(gridE)) ) 
  spec[1,:]=gridE
  for itime in range(ntimes):
   spec[0,:]=times[itime]*1.
   if l==3:
    spec[3:,:]=np.sum(probs[itime,:,:,:],axis=0)
   else:
    spec[3:,:]=probs[itime,l,:,:]
   spec[2,:]=np.sum(spec[3:,:],axis=0)
   np.savetxt(f,np.transpose(spec))
   f.write("\n")

