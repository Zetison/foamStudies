#!/bin/bash
MODEL_NAME="cylinder"
D=1.0
b=64.0
Re=100
t=4.0
NP=2
NT=1
scriptLoc=$(pwd)
outputFolder=$HOME/results/IFEM/cylinder
mkdir -p $outputFolder
cp Cylinder2D_chorin-template.xinp Cylinder2D-template.xinp generate_blockMeshDict_file.c parameters exportResults.py $outputFolder
cd $outputFolder
NAVIERSTOKES=$HOME/kode/IFEM/Apps/IFEM-NavierStokes/r-mpi/bin/NavierStokes
#NAVIERSTOKES=/home/akva/kode/IFEM/Apps/IFEM-NavierStokes/r-mpi/bin/NavierStokes
GENERATOR=$HOME/kode/meshscripts/cylinder/cylinder.py
FORMS="mixed mixed-full subgrid chorin"
FORMS="mixed-full"
GENERATE=$1
POSTPROCESS="-vtf 1 -hdf5"
POSTPROCESS=""
RUN=$2
p_arr=$(seq 1 3)
p_arr=1
MESH_ARR=$(seq 6 7)
MESH_ARR=6

SED_VARIABLES="BDD DIAM LENGTH NRE_0 RE MESH OMEGA_ROT \
							 NP START_TIME END_TIME DELTA_T U_INF AREF LREF \
							 EPSILON_0 K_0 NUT_0 OMEGA_0 NU MODEL_NAME WRITE_INTERVAL PURGE_WRITE STRIDE \
							 V_CLASS P_CLASS V_TYPE V_PC V_PACKAGE P_TYPE P_PC P_PACKAGE \
							 FORM NGAUSS"
							 
if [[ $GENERATE == 1 ]]
then
	for MESH in $MESH_ARR
	do
		nel=$((2**($MESH-1)))
		meshdir=MESH$MESH
		mkdir -p $meshdir
		source ./parameters
		cp generate_blockMeshDict_file.c ./$meshdir/
		pushd $meshdir > /dev/null
		DICTS="generate_blockMeshDict_file.c"
		source $HOME/kode/bashScripts/sedFiles
		g++ generate_blockMeshDict_file.c -o generate_blockMeshDict_file 
		beta=$(./generate_blockMeshDict_file "ifem")
		popd > /dev/null

	  for p in $p_arr
	  do
	    dir=MESH$MESH/$p
	    mkdir -p $dir
	    pushd $dir > /dev/null
	    echo "Generating model for p=$p and MESH=$MESH"
			python3 $GENERATOR --diam=$D --width=$b --front=$b --thickness=$t --back=192.0 --side=$b --height=0.0 --Re=$Re --grad=$beta --nel-bndl=$NRE --nel-circ=$nel --nel-height=1 --order=$((p+1)) --no-outer-graded --out cyl2D
	    popd > /dev/null
	
	
	    for FORM in $FORMS
	    do
				caseFolder=$dir/$FORM
	      mkdir -p $caseFolder
				SOLVER='iterative'
				V_CLASS=petsc
				P_CLASS=$V_CLASS
	      if test "$FORM" = "chorin"
	      then
					V_TYPE=gmres
					V_PC=asm
					V_PACKAGE=petsc
					P_TYPE=$V_TYPE
					P_PC=gamg
					P_PACKAGE=$V_PACKAGE
					NGAUSS=$((p+1))
	        cp Cylinder2D_chorin-template.xinp $caseFolder/Cyl2D.xinp
	      else
					V_TYPE=preonly
					V_PC=lu
					V_PACKAGE=mumps
					P_TYPE=$V_TYPE
					P_PC=$V_TYPE
					P_PACKAGE=$V_PACKAGE
					NGAUSS=$((p+2))
	        cp Cylinder2D-template.xinp $caseFolder/Cyl2D.xinp
	      fi
	    	pushd $caseFolder > /dev/null
				DICTS="Cyl2D.xinp"
				source $HOME/kode/bashScripts/sedFiles
	    	popd > /dev/null
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
  	    pushd MESH$MESH/$p/$FORM > /dev/null
  	    OMP_NUM_THREADS=$NT mpirun -np $NP $NAVIERSTOKES Cyl2D.xinp $POSTPROCESS -petsc -msgLevel 1 | tee Cyl2D.log
  	    popd > /dev/null
				python3 exportResults.py --inputname="MESH"$MESH"/$p/$FORM/Cyl2D_force.dat" --outputname="results_IFEM_"$FORM"_p"$p"_MESH"$MESH".txt" --ifem
  	  done
  	done
	done
fi










