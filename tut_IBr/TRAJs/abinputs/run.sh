if [ ! -d TEMP ]; then
 mkdir TEMP
fi
cd TEMP
cp ../geom0.xyz geom.xyz
echo "3 -1" | model_IBr.exe > model.out

echo "1 0 0
      0 1 0
      0 0 1" > overlap.dat
