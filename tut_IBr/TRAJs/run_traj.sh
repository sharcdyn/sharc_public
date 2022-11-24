#!/bin/bash

acttraj=$SLURM_ARRAY_TASK_ID
if [ $1 ]; then
 acttraj=$1
fi

if [ ! -d TRAJ$acttraj ]; then
 mkdir TRAJ$acttraj
fi
cd TRAJ$acttraj

ln -sf ../abinputs/run.sh .
ln -sf ../abinputs/take.sh .
ln -sf ../abinputs/run_grad.sh .

if [ ! -e in ]; then
 cp ../abinputs/in .
fi

if [ ! -e geom ]; then
 if [ -e ../geoms ]; then
  fichout="../geoms"
 elif [ -e ../geoms.init ]; then
  fichout="../geoms.init"
 fi
 awk -v traj=$acttraj '
  /Traj/{if ($1==traj) {pr=1; nl=-1;pot=$6}}
  {if (NR==1) {nat=$1}
   if (NR>1 && NR<=nat+1) {at[NR-1]=$1;noat[NR-1]=$2;mass[NR-1]=$3}
   if (pr==1) {
    nl+=1
    if (nl>0 && nl<=nat) {print at[nl],noat[nl],$1,$2,$3,mass[nl]}
    if (nl>nat){print}
    if (nl==2*nat) {pr=0}
   }
  }
  END{if (pot*1 !=0 ){print pot}}
 ' $fichout > geom
fi 

sharc.exe < in > out
