#! /bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1:22
#$ -o /data/scratch/hmy117



cd /data/Wolfson-UKBB-Dobson/SMR


./smr_Linux --bfile ../1kg_reference/filtered_chr${SGE_TASK_ID} \
--gwas-summary ./gwas/ms_gwas.ma \
--beqtl-summary ./datasets/westra_eqtl_hg19 \
--out ./results/ms_smr_chr${SGE_TASK_ID} \
--diff-freq-prop 1 \
--thread-num $NSLOTS 

cat ./results/ms_smr_chr${SGE_TASK_ID}\.smr >> ./results/overall_ms_westra_smr_results.txt
rm ./results/ms_smr_chr${SGE_TASK_ID}\.smr*

