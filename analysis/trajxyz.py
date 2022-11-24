def readtraj(fich):
 traj={}
 with open(fich,"r") as f:
  sat=[]
  title=[]
  pos=[]
  vel=[]
  nat=None
  n=-1
  nstep=0
  for line in f:
   kk=line.split()
   if (nat==None): nat=int(kk[0])
   if n==-1 and int(kk[0])!=nat:
    print("Number of atoms cannot be different") 
   if n==0:
    title.append(kk)
    pos0=[]
    vel0=[]
   if n>=1:
    if nstep==0: sat.append(kk[0])
    pos0.append([float(kk[1])/.5292,float(kk[2])/.5292,float(kk[3])/.5292])
    if len(kk)==7:
     vel0.append([float(kk[4])/1000,float(kk[5])/1000,float(kk[6])/1000])
    else:
     vel0.append([ 0.0, 0.0, 0.0 ])
   if n==nat:
    nstep+=1
    n=-1
    pos.append(pos0)
    vel.append(vel0)
   else:
    n+=1
 traj['sat']=sat
 traj['nat']=nat
 traj['title']=title
 traj['pos']=pos
 traj['vel']=vel

 return traj
