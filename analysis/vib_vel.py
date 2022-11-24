import numpy as np
import trajxyz
import sys,os

trajs=[]
dir_trajs=[]
if os.path.isfile("traj.xyz"):
 dir_trajs.append("./")
 trajs.append(trajxyz.readtraj("traj.xyz"))
elif len(sys.argv)>1:
 for i in sys.argv[1:]:
  print("Reading TRAJ ",int(i))
  trajs.append(trajxyz.readtraj("TRAJ"+str(int(i))+"/traj.xyz"))
  dir_trajs.append( "TRAJ"+str(int(i))+"/" )
elif os.path.isfile("valids.dat"):
 val=np.loadtxt("valids.dat")
 for i in val:
  print("Reading TRAJ ",int(i))
  trajs.append(trajxyz.readtraj("TRAJ"+str(int(i))+"/traj.xyz"))
  dir_trajs.append( "TRAJ"+str(int(i))+"/" )
else:
 print("I cannot read the xyz files")
 quit()

print("Number of trajs in memory :",len(trajs))

def velft(traj):
 vel=np.asarray(traj['vel'][:])
 times=[]
 for i in traj['title']:
  times.append(float(i[2]))
 vel0=vel[0,:,:]
 nvel0=np.dot(vel0.ravel(),vel0.ravel())
 nsteps=vel.shape[0]
 nat=vel.shape[1]
 times=np.asarray(times)
 
 corr=np.zeros( (nsteps) )
 for i in range(nsteps):
  corr[i]=np.dot(vel[i,:,:].ravel(),vel0.ravel())
 
 ## FT
 if nvel0>1e-10: corr=corr/nvel0
 corr2=np.fft.fft(corr,norm="ortho")

 dt=times[-1]-times[0]
 dt=dt*41.32/float(times.shape[0]-1)
 dw=2.0*np.pi/dt
 dw=dw/float(nsteps)
 nw=int(nsteps/2)
 corr2=corr2[:nw]
 corr_prob=np.zeros( (nw,2) )
 corr_prob[:,1]=np.abs(corr2)**2
 for i in range(nw):
  corr_prob[i,0]=float(i)*dw
 corr0=np.zeros( (nsteps,2) )
 corr0[:,0]=times
 corr0[:,1]=corr

 return corr_prob,corr2,corr0

specs=[]
nsteps=None
for traj,dire in zip(trajs,dir_trajs):
 spec_prob,spec,corr=velft(traj)
 specs.append(spec_prob)
 np.savetxt(dire+"/ft.dat.gz",spec_prob)
 np.savetxt(dire+"/corr.dat.gz",corr)
 print("Done FT and saved in ",dire+"/ft.dat.gz and ",dire+"/corr.dat.gz")
 if nsteps==None:
  nsteps=spec_prob.shape[0]
 if nsteps>0 and spec_prob.shape[0]!=nsteps:
  nsteps=-1
  print("The trajectory do not have the same number of steps")

if nsteps>0:
 tot_spec=np.zeros( (nsteps,2) )
 for i in specs:
  tot_spec[:,0]=i[:,0]
  tot_spec[:,1]+=i[:,1]
else:
 steps=[]
 for i in specs:
  steps.append(i.shape[0])
 istep=np.argmin(steps)
 nsteps=steps[istep]
 print("Taking the lowest number of steps ",nsteps,istep)
 tot_spec=np.zeros( (nsteps,2) )
 tot_spec[:,0]=specs[istep][:,0]
 for i in specs:
  for j in range(nsteps):
   istep=np.argmin((tot_spec[j,0]-i[:,0])**2)
   tot_spec[j,1]+=i[istep,1]
tot_spec[:,1]=tot_spec[:,1]/float(len(specs))
print("Saving average spectrum in avft.dat.gz ",tot_spec.shape)
np.savetxt("avft.dat.gz",tot_spec)
