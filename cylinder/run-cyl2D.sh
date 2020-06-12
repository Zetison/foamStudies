#!/bin/bash
D=1.0
b=64.0
Re=100
t=4.0
NP=32
NT=1
scriptLoc=$(pwd)
outputFolder=$HOME/results/IFEM/cylinder
mkdir -p $outputFolder
cp Cylinder2D_chorin-template.xinp Cylinder2D-template.xinp generate_blockMeshDict_file.c sed_parameters exportResults.py $outputFolder
cd $outputFolder
#NAVIERSTOKES=$HOME/kode/IFEM/Apps/IFEM-NavierStokes/r-mpi/bin/NavierStokes
NAVIERSTOKES=/home/akva/kode/IFEM/Apps/IFEM-NavierStokes/r-mpi/bin/NavierStokes
GENERATOR=$HOME/kode/meshscripts/cylinder/cylinder.py
FORMS="chorin mixed mixed-full subgrid"
#FORMS="mixed-full"
GENERATE=$1
POSTPROCESS="-vtf 1 -hdf5"
#POSTPROCESS=""
RUN=$2
p_arr=$(seq 1 3)
#p_arr=1
MESH_ARR=$(seq 6 7)

if [[ $GENERATE == 1 ]]
then
	for MESH in $MESH_ARR
	do
		nel=$((2**($MESH-1)))
		nRe=$(echo "5*2^($MESH-6)" | bc -l) 				  # minimum number of elements withing D/sqrt(Re) outside cylinder
		meshdir=MESH$MESH
		mkdir -p $meshdir
		cp generate_blockMeshDict_file.c Cylinder2D_chorin-template.xinp Cylinder2D-template.xinp ./$meshdir/
		pushd $meshdir > /dev/null
		../sed_parameters "cylinder" $NP "ifem" $MESH
		g++ generate_blockMeshDict_file.c -o generate_blockMeshDict_file 
		beta=$(./generate_blockMeshDict_file "ifem")
		popd > /dev/null

	  for p in $p_arr
	  do
	    dir=MESH$MESH/$p
	    mkdir -p $dir
	    pushd $dir > /dev/null
	    echo "Generating model for p=$p and MESH=$MESH"
			python3 $GENERATOR --diam=$D --width=$b --front=$b --thickness=$t --back=192.0 --side=$b --height=0.0 --Re=$Re --grad=$beta --nel-bndl=$nRe --nel-circ=$nel --nel-height=1 --order=$((p+1)) --no-outer-graded --out cyl2D
	    popd > /dev/null
	
	
	    for FORM in $FORMS
	    do
	      mkdir -p $dir/$FORM
	      if test "$FORM" = "chorin"
	      then
					NGAUSS=$((p+1))
	        sed -e "s/SED_FORM/$FORM/g" -e "s/SED_NGAUSS/$NGAUSS/g" $meshdir/Cylinder2D_chorin-template.xinp > $dir/$FORM/Cyl2D-Re$Re.xinp
	      else
					NGAUSS=$((p+2))
	        sed -e "s/SED_FORM/$FORM/g" -e "s/SED_NGAUSS/$NGAUSS/g" $meshdir/Cylinder2D-template.xinp > $dir/$FORM/Cyl2D-Re$Re.xinp
	      fi
	      cp $dir/cyl2D.xinp $dir/cyl2D.g2 $dir/$FORM
	    done
	    rm $dir/cyl2D.xinp $dir/cyl2D.g2
	  done
	done
fi

if [[ $RUN == 1 ]]
then
	for MESH in $MESH_ARR
	do
  	for p in $p_arr
  	do
  	  for FORM in $FORMS
  	  do
  	    pushd Re$Re/$p/$FORM
  	    OMP_NUM_THREADS=$NT mpirun -np $NP $NAVIERSTOKES Cyl2D-Re$Re.xinp $POSTPROCESS -petsc -msgLevel 1 | tee Cyl2D-Re$Re.log
  	    popd
				python3 exportResults.py --inputname="MESH"$MESH"/$p/$FORM/Cyl2D_force.dat" --outputname="results_IFEM_"$FORM"_p"$p"_MESH"$MESH".txt" --ifem
  	  done
  	done
	done
fi










