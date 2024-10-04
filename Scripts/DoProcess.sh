#!/bin/bash

## PC
export MPI="mpirun -n 16"
export OMP_NUM_THREADS=16

## Case
export KaF=0


## runtime vars
set -x
export i_H2=90
export PROG_H2=`echo "0.01*${i_H2}" | bc -l`

if [ ${KaF} -eq 0 ] || [ ${KaF} -eq 1 ]
then
    export NPTS=100000
elif [ ${KaF} -eq 4 ]
then
    export NPTS=150000
elif [ ${KaF} -eq 12 ]
then
    export NPTS=200000
elif [ ${KaF} -eq 36 ]
then
    export NPTS=250000
fi
export NCELLSPERLF=16
export JHI=`echo "12*${NCELLSPERLF}" | bc -l`
export NRK=`echo "2*${JHI} + 1" | bc -l`
export HRK="0.25"
export GROW=`echo "3*${NCELLSPERLF} + 2" | bc -l`


for PLT in plt00000_FLOAT; do
#    ${MPI} ./plotProg3d.gnu.MPI.ex inputs.process infile= ${PLT} vars="temp Y(H2)" normalisePltName= ${PLT};
    #${MPI} ./plotYtoPhi3d.gnu.MPI.ex inputs.process infile= ${PLT};
#    ${MPI} ./grad3d.gnu.MPI.ex inputs.process infile= ${PLT} gradVar="Y(H2)" normalisePltName= ${PLT};
#    ${MPI} ./combinePlts3d.gnu.MPI.ex inputs.process infiles= ${PLT} ${PLT}_prog vars=x_velocity y_velocity z_velocity "Y(H2)_prog" outfile= ${PLT}_vel_prog;
#    ${MPI} ./curvature3d.gnu.MPI.ex inputs.curvature progressName= "Y(H2)_prog" infile= ${PLT}_vel_prog outfile= ${PLT}_vel_prog_K
#    ${MPI} ./combinePlts3d.gnu.MPI.ex inputs.process infiles= ${PLT} ${PLT}_vel_prog_K vars=x_velocity y_velocity z_velocity temp HeatRelease H2_ConsumptionRate "Y(H2)" "Y(H2)_prog" "MeanCurvature_Y(H2)_prog" "GaussianCurvature_Y(H2)_prog" "StrainRate_Y(H2)_prog" outfile= ${PLT}_combine;
#    ${MPI} isosurface3d.gnu.MPI.ex inputs.process \
#	   comps=7 \
#	   isoCompName= "Y(H2)_prog" \
#	   isoVal=${PROG_H2} \
#	   infile= ${PLT}_combine \
#	   outfile_base=${PLT}_iso
#    ./decimateMEF3d.gnu.ex -t ${NPTS} ${PLT}_iso.mef > ${PLT}_iso_qslim.mef;


    ${MPI} partStream3d.gnu.MPI.ex infile= ${PLT}_combine \
	   isofile= ${PLT}_iso_qslim.mef \
	   writeStreams=0 writeParticles=0 \
	   hRK= ${HRK} Nsteps=${NRK} \
	   vars= x_velocity y_velocity z_velocity "FlameNormalX_Y(H2)_prog" "FlameNormalY_Y(H2)_prog" "FlameNormalZ_Y(H2)_prog"
	   outfile= ${PLT}
##infile= pltX vectorField= gradX_x gradX_y (gradX_z) vars= field1 field2 ... (cSpace= 1, def= 0) (hRK= 0.4, def= 0.1) (Nsteps= 100, def= 50)    


    
##	srun ./plotPhi3d.gnu.x86-rome.MPI.ex inputs.process infile= ../${PLT} outfile=${PLT}_phi
##	mv ${PLT}_prog ${PLT}_phi
##	sleep 2
##	srun ./plotProg3d.gnu.x86-rome.MPI.ex inputs.process infile= ../${PLT}
##	#srun ./plotQoverRhoCP3d.gnu.x86-rome.MPI.ex inputs.process infile= ../${PLT} outfile=${PLT}_QRhoCP 
##	srun ./combinePlts3d.gnu.x86-rome.MPI.ex inputs.process infiles=../${PLT} ${PLT}_prog outfile=${PLT}_prog_comb vars=x_velocity y_velocity z_velocity temp I_R\(H2\) HeatRelease Y\(H2\) prog_H2
##        srun ./curvature3d.gnu.x86-rome.MPI.ex inputs.curvature progressName= prog_H2 infile= ${PLT}_prog_comb outfile= ${PLT}_curv
##	srun ./combinePlts3d.gnu.x86-rome.MPI.ex inputs.process infiles=${PLT}_prog_comb ${PLT}_curv ${PLT}_gt ${PLT}_phi ${PLT}_QoverRhoCP outfile=${PLT}_curv_comb vars=x_velocity y_velocity z_velocity temp I_R\(H2\) HeatRelease Y\(H2\) prog_H2 MeanCurvature_prog_H2 GaussianCurvature_prog_H2 StrainRate_prog_H2 FlameNormalX_prog_H2 FlameNormalY_prog_H2 FlameNormalZ_prog_H2 phi QoverRhoCP VelFlameNormal \|\|gradtemp\|\|
##	srun ./isosurface3d.gnu.x86-rome.MPI.ex inputs.process comps= 7 isoCompName= prog_H2 isoVal= ${PROG_H2} infile=${PLT}_prog_comb outfile= ${PLT}_prog_comb_prog_H2_${i_H2}.mef
##	srun ./surfMEFtoDAT3d.gnu.x86-rome.MPI.ex inputs.process infile= ${PLT}_prog_comb_prog_H2_0.9.mef outfile= ${PLT}_prog_comb_prog_H2_${PROG_H2}.dat 
	#./qslim3d.gnu.ex -t ${NPTS} ${PLT}_prog_H2_${i_H2}.mef > ${PLT}_prog_H2_${i_H2}_${NPTS}.mef
##	./decimateMEF3d.gnu.x86-rome.ex -t $NPTS ${PLT}_prog_comb_prog_H2_0.9.mef > ${PLT}_prog_comb_prog_H2_${PROG_H2}.$NPTS.mef	
##	srun ./partStream3d.gnu.x86-rome.MPI.ex infile= ${PLT}_curv_comb isofile= ${PLT}_prog_comb_prog_H2_${PROG_H2}_${NPTS}.mef \
##		writeStreams=1 writeParticles=1 is_per= 1 1 0 outfile= ${PLT}_prog_H2_${i_H2}_${NPTS} Nsteps= ${JHI} 
	#srun ./
done

#for j in $basename/Ka*/chi*; do
#dirname="$j"
#cp $TOP/convertPlt2Float3d.gnu.x86-rome.FLOAT.MPI.OMP.ex $dirname
#echo $dirname
#cd $dirname
#for PLT in plt????0; do 
#
#    srun ./convertPlt2Float3d.gnu.x86-rome.FLOAT.MPI.OMP.ex infile=${PLT} outfile=${PLT}f \
#	 remove= RhoRT I_R\(H\) I_R\(O\) I_R\(O2\) I_R\(OH\) I_R\(H2O\) I_R\(HO2\) \
#	 I_R\(CH2\) I_R\(CH2\(S\)\) I_R\(CH3\) I_R\(CO\) I_R\(CO2\) I_R\(HCO\) \
#	 I_R\(CH2O\) I_R\(CH3O\) I_R\(C2H4\) I_R\(C2H5\) I_R\(C2H6\) I_R\(N2\) I_R\(AR\) \
#	 FunctCall \
#	 >> output.reduce.$SLURM_JOBID 2>>output.reduce.$SLURM_JOBID
#
#done
#done
