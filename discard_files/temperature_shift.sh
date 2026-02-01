#!/bin/bash
#PBS -q c56m256g
#PBS -l mem=8gb,nodes=1:ppn=2,walltime=2000:00:00
#PBS -e temperature_shift.e
#PBS -o temperature_shift.o
#HSCHED -s hschedd

SCRIPT=/p300s/zhaowm_group/tangbx/idog/test/script

cd /gpfs/zhaowm_group/tangbx/cell/DEG

source activate /p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4
python $SCRIPT/mergeRSEM.py -i /gpfs/zhaowm_group/tangbx/cell/projects_results/rnseq -l /gpfs/zhaowm_group/tangbx/cell/DEG/temperature_shift/temperature_shift.list -o /gpfs/zhaowm_group/tangbx/cell/DEG/temperature_shift/temperature_shift.count.tsv -c expected_count
conda deactivate
