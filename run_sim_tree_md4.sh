#!/bin/bash
#  run_sim_tree_md4.sh
#
#  Created by A Lyne on 03/03/2017.
#
#PBS -l nodes=1:ppn=1

# set max wallclock time
#PBS -l walltime=40:00:00

# set memory requirements
#PBS -l mem=3gb

# set name of job
#PBS -N run_repeats_md4

# get output as job is running
#PBS -k oe

# number of tasks to submit
#PBS -t 1-60

# start job from the directory it was submitted
cd $PBS_O_WORKDIR

R_EXEC="/bioinfo/local/build/Centos/R/R-3.4.0/bin/Rscript"
R_FILES1="run_repeats_inference_cluster.R"
R_FILES2="run_post_paup.R"
R_FILES3="run_post_revbayes.R"
R_FILES4="record_NAs.R"
TREE=5
TUNE=2000
BURN=2000
ITER=8000
MDSCEN=4 #missing data scenario 1=no missing data, 2=uniform [0,1] missing data, 3=10% missing data, 4=20% missing data, 5=30% missing data, 6=scrnaseq scenario
mcmc_trees1=${PBS_O_WORKDIR}/mcmc_trees_"$PBS_ARRAYID"_FK_md_"$MDSCEN".trees
rm nst_"$PBS_ARRAYID"_md_"$MDSCEN".txt

for ((i=1; i<=2; i++)); do
    for ((j=1; j<=5; j++)); do
        echo "i=$i j=$j"
        #run first r file to do distance based reconstruction
        $R_EXEC $R_FILES1 $PBS_ARRAYID $TREE $i $j $MDSCEN $PBS_O_WORKDIR
        NST=$(awk 'NR==1' nst_"$PBS_ARRAYID"_md_"$MDSCEN".txt)
        echo $NST
        if [ "$NST" != 0 ]; then
            #change relevant parts of nexus file for PAUP
            #replace PROTEIN with STANDARD and add lines at end of file
            sed -i 's/PROTEIN/STANDARD SYMBOLS="A~Z"/g' ms_sim_data_paup_"$PBS_ARRAYID"_md_"$MDSCEN".nex
            echo -e "BEGIN ASSUMPTIONS;\nOPTIONS DEFTYPE=ORD;\nENDBLOCK;\n\nBEGIN PAUP;\nset autoclose=yes;\nset criterion=parsimony;\nset storebrlens=yes;\nset increase=auto;\ncondense collapse=no;\nhsearch addseq=random nreps=5 swap=tbr hold=1 rearrlimit=1000000;\nsavetrees file=pars_paup_output${PBS_ARRAYID}_md_${MDSCEN}.nex format=altnex brlens=yes from=1 to=1;\nENDBLOCK;" >> ms_sim_data_paup_"$PBS_ARRAYID"_md_"$MDSCEN".nex
            #run paup and execute file
            STARTTIME=$(date +%s)
            ./paup4a166_centos64 -n ms_sim_data_paup_"$PBS_ARRAYID"_md_"$MDSCEN".nex
            ENDTIME=$(date +%s)
            echo $(($ENDTIME - $STARTTIME)) > timing_"${PBS_ARRAYID}"_md_"${MDSCEN}".txt
            
            #run second R file to compute tree distances
            $R_EXEC $R_FILES2 $PBS_ARRAYID $TREE $i $j $MDSCEN $PBS_O_WORKDIR
            rm timing_"${PBS_ARRAYID}"_md_"${MDSCEN}".txt

            #run revbayes 1
            rb_command="int_seed = $PBS_ARRAYID;
            n_state=$NST;
            data_file=\"${PBS_O_WORKDIR}/ms_sim_data_rb_${PBS_ARRAYID}_md_${MDSCEN}.tsv\";
            tree_file=\"${PBS_O_WORKDIR}/bayes_tree_${PBS_ARRAYID}_md_${MDSCEN}.tre\";
            mcmc_trees=\"$mcmc_trees1\";
            mcmc_tune=$TUNE;
            mcmc_burn=$BURN;
            mcmc_iter=$ITER;
            source(\"MS_bayes_freeK_discrete.Rev\");"
            STARTTIME=$(date +%s)
            echo $rb_command | rb
            ENDTIME=$(date +%s)
            echo $(($ENDTIME - $STARTTIME)) > timing_"${PBS_ARRAYID}"_md_"${MDSCEN}".txt

            #run R file to compute distances
            $R_EXEC $R_FILES3 $PBS_ARRAYID $TREE $i $j $mcmc_trees1 $MDSCEN $PBS_O_WORKDIR
            rm timing_"${PBS_ARRAYID}"_md_"${MDSCEN}".txt

            #delete files
            rm ms_sim_data_paup_"$PBS_ARRAYID"_md_"$MDSCEN".nex ms_sim_data_rb_"$PBS_ARRAYID"_md_"$MDSCEN".tsv pars_paup_output"$PBS_ARRAYID"_md_"$MDSCEN".nex map_tree_"$PBS_ARRAYID"_md_"$MDSCEN".txt ml_tree_"$PBS_ARRAYID"_md_"$MDSCEN".txt nst_"$PBS_ARRAYID"_md_"$MDSCEN".txt nj_tree*_"$PBS_ARRAYID"_"$i"_"$j"_"$MDSCEN".new bal_tree*_"$PBS_ARRAYID"_"$i"_"$j"_"$MDSCEN".new ml_tree_"$PBS_ARRAYID"_md_"$MDSCEN".new map_tree_"$PBS_ARRAYID"_md_"$MDSCEN".new pars_paup_output"$PBS_ARRAYID"_md_"$MDSCEN".new mcmc_trees_"$PBS_ARRAYID"_FK_md_"$MDSCEN".trees pars_paup_output"$PBS_ARRAYID"_md_"$MDSCEN".nex balcos_tree_"$PBS_ARRAYID"_"$i"_"$j"_"$MDSCEN".new njcos_tree_"$PBS_ARRAYID"_"$i"_"$j"_"$MDSCEN".new
        else
            #record NA in results matrix
            $R_EXEC $R_FILES4 $PBS_ARRAYID $TREE $i $j $MDSCEN
            rm nst_"$PBS_ARRAYID"_md_"$MDSCEN".txt
        fi
    done
done



