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
--beqtl-summary /data/Wolfson-UKBB-Dobson/SMR/datasets/GTEx_V7_cis_eqtl_summary_lite/Whole_Blood_1e-05 \
--out ms_smr_chr${SGE_TASK_ID} \
--diff-freq-prop 1 \
--thread-num $NSLOTS

cat ms_smr_chr${SGE_TASK_ID}\.smr >> overall_ms_smr_gtex_results.txt
rm ms_smr_chr${SGE_TASK_ID}\.smr*



