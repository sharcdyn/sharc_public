lista=` ls -1 -d  TRAJ* ` 
echo $lista
for i in $lista; do
 if [ -d $i ]; then
  cd $i
  sharc_internal.exe < ../internal.in
  cd ..
 fi
done
