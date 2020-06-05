#!/bin/bash
module load openmpi-x86_64
scl enable devtoolset-8 bash
D=1.0
b=64.0
Re=100
MESH=7
nel=$((2**($MESH-1)))
nRe=10
t=4.0
NP=2
NT=2
STRIDE=10
scriptLoc=$(pwd)
outputFolder=$HOME/hugeFiles/IFEM/cylinder
cp Cylinder2D_chorin-template.xinp Cylinder2D-template.xinp generate_blockMeshDict_file.c sed_parameters exportResults.py $outputFolder
cd $outputFolder

./sed_parameters "cylinder" $NP "ifem" $MESH
g++ generate_blockMeshDict_file.c -o generate_blockMeshDict_file 
beta=$(./generate_blockMeshDict_file "ifem")
#NAVIERSTOKES=$HOME/kode/IFEM/Apps/IFEM-NavierStokes/r-mpi/bin/NavierStokes
NAVIERSTOKES=/home/akva/kode/IFEM/Apps/IFEM-NavierStokes/r-mpi/bin/NavierStokes
GENERATOR=$HOME/kode/meshscripts/cylinder/cylinder.py
#FORMS="chorin mixed mixed-full subgrid"
FORMS="mixed-full"
GENERATE=$1
RUN=$2

if [[ $GENERATE == 1 ]]
then
  for p in `seq 1 1`
  do
    dir=Re$Re/$p
    mkdir -p $dir
    pushd $dir > /dev/null
    echo "Generating model for p=$p"
		python3 $GENERATOR --diam=$D --width=$b --front=$b --thickness=$t --back=192.0 --side=$b --height=0.0 --Re=$Re --grad=$beta --nel-bndl=$nRe --nel-circ=$nel --nel-height=1 --order=$((p+1)) --no-outer-graded --out cyl2D
    popd > /dev/null


    NGAUSS=3
    if test $p -ge 2
    then
      NGAUSS=4
    fi
    for FORM in $FORMS
    do
      mkdir -p $dir/$FORM
      if test "$FORM" = "chorin"
      then
        sed -e "s/SED_FORM/$FORM/g" -e "s/SED_NGAUSS/$NGAUSS/g" Cylinder2D_chorin-template.xinp > $dir/$FORM/Cyl2D-Re$Re.xinp
      else
        sed -e "s/SED_FORM/$FORM/g" -e "s/SED_NGAUSS/$NGAUSS/g" Cylinder2D-template.xinp > $dir/$FORM/Cyl2D-Re$Re.xinp
      fi
      cp $dir/cyl2D.xinp $dir/cyl2D.g2 $dir/$FORM
    done
    rm $dir/cyl2D.xinp $dir/cyl2D.g2
  done
fi

if [[ $RUN == 1 ]]
then
  for p in `seq 1 1`
  do
    for FORM in $FORMS
    do
      pushd Re$Re/$p/$FORM
      OMP_NUM_THREADS=$NT mpirun -np $NP $NAVIERSTOKES Cyl2D-Re$Re.xinp -vtf 1 -hdf5 -petsc -msgLevel 1 | tee Cyl2D-Re$Re.log
      popd
			python3 exportResults.py --inputname="Re"$Re"/$p/$FORM/Cyl2D_force.dat" --outputname="results_IFEM_"$FORM"_p"$p".txt"
    done
  done
fi










