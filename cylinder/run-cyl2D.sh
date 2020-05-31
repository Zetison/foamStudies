#!/bin/bash
D=1.0
b=64.0
Re=100
M=7
N=$((2**($M-1)))
N2=$((2*N))
end_time=1000    # End time for simulation
delta_t=0.05   # Discretization time step
nRe=10
t=4.0
scriptLoc=$(pwd)
outputFolder=$HOME/hugeFiles/IFEM/cylinder
cp Cylinder2D_chorin-template.xinp $outputFolder
cp Cylinder2D-template.xinp $outputFolder
cd $outputFolder 
beta=$(grep -oP '(?<=beta )[0-9]+\.[0-9]+' ~/OpenFOAM/OpenFOAM-v1912/run/cylinder/constant/parameters)
NAVIERSTOKES=$HOME/kode/IFEM/Apps/IFEM-NavierStokes/r-mpi/bin/NavierStokes
GENERATOR=$HOME/kode/meshscripts/cylinder/cylinder.py
NP=4
NT=2

STRIDE=10

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
		python3.7 $GENERATOR --diam=$D --width=$b --front=$b --thickness=$t --back=192.0 --side=$b --height=0.0 --Re=$Re --grad=$beta --nel-bndl=$nRe --nel-circ=$N --nel-height=1 --order=$((p+1)) --no-outer-graded --out cyl2D
    popd > /dev/null


    NGAUSS=3
    if test $p -ge 2
    then
      NGAUSS=4
    fi
    for FORM in chorin # mixed mixed-full subgrid
    do
      mkdir -p $dir/$FORM
      if test "$FORM" = "chorin"
      then
        sed -e "s/DELTA_T/$delta_t/g" -e "s/END_TIME/$end_time/g" -e "s/STRIDE/$STRIDE/g" -e "s/FORM/$FORM/g" -e "s/REYNOLDS/$Re/g" -e "s/NGAUSS/$NGAUSS/g" Cylinder2D_chorin-template.xinp > $dir/$FORM/Cyl2D-Re$Re.xinp
      else
        sed -e "s/DELTA_T/$delta_t/g" -e "s/END_TIME/$end_time/g" -e "s/STRIDE/$STRIDE/g" -e "s/FORM/$FORM/g" -e "s/REYNOLDS/$Re/g" -e "s/NGAUSS/$NGAUSS/g" Cylinder2D-template.xinp > $dir/$FORM/Cyl2D-Re$Re.xinp
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
    for FORM in chorin # mixed mixed-full subgrid
    do
      pushd Re$Re/$p/$FORM
      OMP_NUM_THREADS=$NT mpirun -np $NP $NAVIERSTOKES Cyl2D-Re$Re.xinp -vtf 1 -hdf5 -petsc -msgLevel 1 | tee Cyl2D-Re$Re.log
      popd
    done
  done
fi










