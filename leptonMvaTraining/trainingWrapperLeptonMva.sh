#!/bin/bash
source ../scripts/setCMSSW.sh

for leptonFlavor in muon electron; do
    for algorithm in BDT MLP DNN; do
        script=train_${leptonFlavor}_${algorithm}.sh
        > $script
        makeSubmit ${script} $PWD
        echo "./trainLeptonMva $leptonFlavor $algorithm" >> $script
        
        if [ "$algorithm" == "DNN" ]; then
            qsub $script -lnodes=1:ppn=16 -l walltime=80:00:00
        else 
            qsub $script -l walltime=40:00:00
        fi
    done
done    
